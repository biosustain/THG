#!/usr/bin/python
import copy
import logging
import os
import pickle
import re
import urllib
from functools import reduce
from typing import Dict, List

import cobra
from cobra.io import read_sbml_model
import dill

# reimports for type hints
from class_generate_database import *
from class_generate_database import compound as CompoundType
from class_generate_database import gene as GeneType
from class_generate_database import reaction as ReactionType
from equations_generate_database import * 
from functions_generate_database import *
from pattern_generate_database import *


def cobra_reconstruction(
    model_name: str,
    model_id: str,
    metabolite_list: Dict[List, CompoundType],
    reaction_list: Dict[List, ReactionType],
    gene_list: Dict[List, GeneType],
    pathways: Dict[str, str],
    location_dict: Dict[str, str],
    metabolite_equivalent: Dict[str, str],
    metabolite_list_general: Dict[List, CompoundType],
) -> cobra.Model:
    """Reconstruction of gathered information given by using cobrapy.

    Parameters
    ----------
    model_name: str
    model_name: id
    metabolite_list: Dict[List, compound]
    reaction_list: Dict[List, reaction]
    gene_list: Dict[List, gene]
    pathways: Dict[str, str]
        Map from Pathway to Reaction identifiers (PathNameRxn).
    location_dict: Dict[str, str]
        generated from the Comparment_Cl. This stores human readable
        compartment names to comparment identifiers. The extracted identifiers
        in the reactions and metabolites contain the id and the comparment
        name, so this mapping is necessary to achieve a proper
        SBML-compatible identifier (no spaces).
    """
    model = cobra.Model(model_id or model_name, model_name or model_id)
    location_dict = {k.lower(): v for k, v in location_dict.items()}
    # metabolites
    compounds = [
        cobra.Metabolite(
            compound.ID2() + "_" + location_dict.get(compound.Subcel.lower()),
            compound.Formula1(),
            compound.Name(),
            float(compound.charge()),
            location_dict.get(compound.Subcel.lower()),
        )
        for iden, compound in metabolite_list.items()
    ]    
    model.add_metabolites(compounds)
    # store mapping (met identifier -> met.id in model) for reaction section
    met_mapping = {}
    # metabolite annotation
    for iden, compound in metabolite_list.items():
        model_met = model.metabolites.get_by_id(
            compound.ID2() + "_" + location_dict.get(compound.Subcel.lower())
        )
        met_mapping[iden] = model_met.id
        annotation = {
            "pubchem.compound": compound.PubChem(),
            "chebi.compound": compound.CheBI(),
            "glycomedb": compound.GlyDB(),
            "jcggdb": compound.JCGGDB(),
            "inchi": compound.inchi(),
            "inchikey": compound.inchikey(),
            "lipidbank": compound.LipidBank(),
            "lipidmaps": compound.LIPIDMAPS(),
        }
        alt_formulas = [compound.Formula2(), compound.Formula3(), compound.Formula4()]
        alt_formula = [
            form for form in alt_formulas if form != model_met.formula and form
        ]
        if alt_formula:
            model_met.annotation["glycan_formula"] = alt_formula[0]

        # only add non-empty annotation
        model_met.annotation = {k: v for k, v in annotation.items() if v}

    for x in model.metabolites: # replace glycan formula by a sbml suitable format
        if x.id[0] == 'G':
            try:
                print(x.id)
                compartment = [y[0] for y in location_dict.items() if y[1] in x.id.split('_')[1]][0]
                xth_metabolite_id = x.id.split('_')[0] + '_' + compartment
                x.formula = metabolite_list[xth_metabolite_id].Formula4()
            except Exception:
                continue

    for x in model.metabolites: # eliminate potential discrepancies between metabolite id and reaction compounds ids
        if x.id.split("_")[0] in metabolite_equivalent.keys():
            x.id = metabolite_list_general[metabolite_equivalent[x.id.split("_")[0]]].ID1() + "_" + x.compartment.lower()

    def normalize_id(reac_id: str):
        iden, comp_desc = reac_id.split("_")
        comp = location_dict[comp_desc.lower()]
        return f"{iden}_{comp}"

    from gpr.ast_gpr import sanitize_gpr
    # reactions
    reactions = [
        cobra.Reaction(
            normalize_id(iden), reac.Name(), "", 0 if reac.Termodyn() else -1000, 1000
        )
        for iden, reac in reaction_list.items()
    ]
    model.add_reactions(reactions)
    for iden, rxn in reaction_list.items():
        print(iden)
        reac = model.reactions.get_by_id(normalize_id(iden))
        rxn_id = rxn.ID
        if "_" in rxn_id:
            kegg_id = rxn_id.split("_")[0]
            comp_id = "_" + location_dict[rxn_id.split("_")[-1].lower()]
        else:
            kegg_id = rxn_id
            comp_id = ""
        try:
            # this may fail after serialization because these are overwritten
            # at runtime via a capturing lambda
            # substrates might come with positive coefficients
            substrates = [
                (-abs(convert_to_float(subs[0])), subs[1], subs[2]) for subs in rxn.Substrate()
            ]
            products = [
                (abs(convert_to_float(prod[0])), prod[1], prod[2]) for prod in rxn.Product()
            ]
            reac_compounds = products + substrates
        except Exception:
            # these metabolites do not have an specified comparment!
            # substrates might come with positive coefficients
            substrates = [(-abs(convert_to_float(subs[0])), subs[1], subs[2]) for subs in rxn.subs]
            products = [(abs(convert_to_float(prod[0])), prod[1], prod[2]) for prod in rxn.prods] 
            reac_compounds = products + substrates
        metabolites = {
            met_mapping[met[2]]
            if met[2] in met_mapping
            else f"{met[2]}{comp_id}": float(met[0])
            for met in reac_compounds
        }
        # there may be some metabolites that were not passed in metabolite list
        new_mets = []
        for met_id in metabolites:
            if met_id not in model.metabolites:
                # first check if we have them in a different comparment
                met_root, met_comp = met_id.split("_")[0], met_id.split("_")[1]
                maybe_mets = [m for m in model.metabolites if m.id.startswith(met_root)]
                if maybe_mets:
                    model_met = maybe_mets[0]
                    LOGGER.warning(
                        f"Reactant '{met_root}' was added in a different comparment '{met_comp}'."
                    )
                    new_met = cobra.Metabolite(
                        met_id,
                        model_met.formula,
                        model_met.name,
                        model_met.charge,
                        # might be a new compartment!
                        met_comp,
                    )
                    new_met.annotation = model_met.annotation
                else:
                    LOGGER.warning(
                        f"Reactant {met_id} was not found in any compartment. Creating new one!"
                    )
                    new_met = cobra.Metabolite(
                        met_id,
                        compartment=met_comp,
                    )
                new_mets.append(new_met)
        if new_mets:
            model.add_metabolites(new_mets)

        reac.add_metabolites(metabolites)
        ec = rxn.EC()
        reac.annotation = {"kegg.reaction": kegg_id, "ec-code": ec[0] if ec else ""}
        sgpr, gpr = rxn.GPR[0], rxn.GPR[1].replace("[", "").replace("]", "")
        if sgpr and sgpr != "[]":
            if "or" in gpr and "and" not in gpr:
                # GPRs scrapped from Kegg are added with ORs and the genes
                # may be repeated so they have to be deduplicated
                # TODO(carrascomj): should come from getGPR / getLocation
                gpr = " or ".join({gene for gene in gpr.split(" or ") if gene})
            reac.gene_reaction_rule = sanitize_gpr(gpr)
            reac.annotation["sGPR"] = sgpr
        reac.id = kegg_id + comp_id
    # add a group per pathway
    model.add_groups([cobra.core.Group(group, group) for group in pathways])
    for group, members in pathways.items():
        # the members are the reactions in each pathway
        # TODO(carrascomj): reaction ids coming from paths are not in compartments
        model.groups.get_by_id(group).add_members(
            reduce(
                lambda x, y: x + y,
                [model.reactions.query(member) for member in members.split()],
                [],
            )
        )
    # gene annotation (genes were added with the GPRs)
    pat_enstp = re.compile("ENS[TP][0-9]+")
    for iden, gene in gene_list.items():
        print(iden)

        if not iden in model.genes:
            LOGGER.warning(f"Gene '{iden}' was not found. Creating new one!")
            model_gene = cobra.Gene(
                    iden,
                    gene.Name().replace("[", "").replace("]", "").replace("-", ""),
                )
            
        else:
            model_gene = model.genes.get_by_id(iden)
            model_gene.name = gene.Name().replace("[", "").replace("]", "").replace("-", "")

        # there may be other ensembl genes, we have to query them
        ensembl = gene.Ensg()
        ensembl_genes = []
        if ensembl:
            try:
                retrieved_ids = str(
                    urllib.request.urlopen(
                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=" + ensembl
                    ).read()
                )  # only add non-empty annotation
                ensembl_genes = [
                    str(x) for x in list(set(pat_enstp.findall(retrieved_ids)))
                ] + [ensembl]
            except Exception:
                ensembl_genes = []
    
        model_gene.annotation = {
            k: v
            for k, v in {
                "ensembl": ensembl_genes,
                "ncbigene": gene.Entrez(),
                "uniprot": gene.Uniprot(),
                "hgcn.symbol": model_gene.name,
            }.items()
            if v
        }
    return model


LOGGER = logging.getLogger(__name__)
session = setup_biocyc_session()

#### Initial Parameters
ListOfPaths = "files/human_kegg_pathways.txt"  # in the current folder
ModelCompounds = "files/extra_compounds.txt"  # in the current folder
ExtraFormula = "files/extra_formula.txt"
ModelReactions = ""
ModelGenes = ""
EnsblDB = "files/ensembl"  # From Ensembl database: ensembl gene ID vs Entrez vs Name.
time = 20  # Time to download url: Parameter defined in function getHtml
Path = (
    open(ListOfPaths, "r").read().split("\n")
)  # From analysis using metaboanalyst.
Compound = (
    open(ModelCompounds, "r").read().split("\n")
)  
EF = [_f for _f in open(ExtraFormula, "r").read().split("\n") if _f]
Output = "model/Human_Database.xml"  # Output model
ModID = Output
ModName = Output
variablesFile = "files/model_variables.pkl"  # File where the working environment is saved
specialCompounds = "files/special_compounds.txt"  # File where we save the IDs of the compounds with a (group)n in their formula
open(specialCompounds, "w").close()  # Erase or create the file


# Dictionary with extra compounds that can be added to mass balance the metabolic reactions
extra_compound = {'H':'C00080','H2O':'C00001','Fe':'C00023','Na':'C01330','Ca':'C00076','K':'C00238','F':'C00023','R':'C00000','X':'C0000X'}


#### Initial List and dictionaries
PathList = {}
RxnList = {}
MetList = {}
GPRList = {}
MetEquiv = {}
PathNameRxn = {}
RxnEquiv = {}
PathIdent = []
RxnIdent = []
MetIdent = []
GPRIdent = []
RxnIDList = []
MetIDList = []


#### List and dictionaries for the subcelular location annotation
CSL_ID = {
"extracellular": "e",
"peroxisome": "x",
"mitochondria": "m",
"cytosol": "c",
"lysosome": "l",
"endoplasmic reticulum": "r",
"golgi apparatus": "g",
"nucleus": "n",
"inner mitochondria": "i",
}  # to keep the consistency between the DB and the initial compartments in Human1
CSL_ID = compartment_file_to_dict()
CSL_ID = dict((k.lower(), v.lower()) for k, v in CSL_ID.items())
listOfID = []
LocVar = {}
RxnList_CL = {}
MetList_CL = {}
RxnIdent_CL = []
MetIdent_CL = []
Compartment_CL = []
RxnList_Subcel = []

######### Pathways ###########
i = 0
while i < len(Path) - 1:

    #### Build network
    PathID = Path[i].split("\t")[0]
    PathName = Path[i].split("\t")[1]
    PathURL = "http://www.kegg.jp/kegg-bin/download?entry=" + PathID + "&format=kgml"
    PathReferer = "https://www.kegg.jp/kegg-bin/show_pathway?" + PathID
    PathList[PathID] = pathway(PathURL, time, PathID, PathReferer, PathName)
    if PathList[PathID].Compounds() or not PathList[PathID].Compounds():
        print(
            PathName
            + ": defined in human"
            + "("
            + str(i + 1)
            + "/"
            + str(len(Path) - 1)
            + ")"
        )
        PathNameRxn[PathName] = ""
        j = 0
        while j < len(PathList[PathID].Reactions()):
            try:
                RxnID = PathList[PathID].Reactions()[j][0][0]
                if (
                    not RxnID in RxnIdent and not RxnID in RxnEquiv
                ):

                    ######### Define New Reaction ###########
                    RxnIdent = RxnIdent + [RxnID]
                    RxnURL = PathList[PathID].Reactions()[j][1]
                    RxnTermDyn = PathList[PathID].Reactions()[j][0][1]
                    RxnList[RxnID] = reaction(RxnURL, time, RxnID, PathName, RxnTermDyn)

                    ######### Check if all the compounds in the jth reaction are in the compound list ###########
                    RxnCmp = [x[2] for x in RxnList[RxnID].Substrate()] + [
                        x[2] for x in RxnList[RxnID].Product()
                    ]
                    c = 0
                    while c < len(RxnCmp):
                        CompID = RxnCmp[c]
                        if not CompID in MetIdent and not CompID in MetEquiv:
                            MetIdent = MetIdent + [CompID]
                            if CompID[0] == "C":
                                CompURL = (
                                    "http://www.kegg.jp/dbget-bin/www_bget?cpd:" + CompID
                                )
                            if CompID[0] == "G":
                                CompURL = (
                                    "http://www.kegg.jp/dbget-bin/www_bget?gl:" + CompID
                                )
                            MetList[CompID] = compound(
                                CompURL, CompID, time, EF, specialCompounds
                            )
                            if MetList[CompID].ID1() != MetList[CompID].ID2():
                                MetIdent[len(MetIdent) - 1] = MetList[CompID].ID1()
                                MetEquiv[CompID] = MetList[CompID].ID1()
                                MetList[MetList[CompID].ID1()] = copy.deepcopy(
                                    MetList[CompID]
                                )  # Change the reference in the dictionary to account for the 1th ID
                                del MetList[CompID]
                        c = c + 1

                    ######### Define Substrates, Products and New Compounds ###########
                    # Evaluate the relation between substrates and products #
                    Rxn = getRxncons(
                        RxnList[RxnID],
                        time,
                        MetEquiv,
                        MetList,
                        MetIdent,
                        EF,
                        specialCompounds,
                    )
                    RxnList[RxnID] = copy.deepcopy(Rxn)
                    
                    # Check reaction ID
                    if RxnID != RxnList[RxnID].ID:
                        RxnEquiv[RxnID] = RxnList[
                            RxnID
                        ].ID  # Glycan Reaction : Compound Reaction
                        RxnIdent[len(RxnIdent) - 1] = RxnList[RxnID].ID
                        RxnList[RxnList[RxnID].ID] = copy.deepcopy(
                            RxnList[RxnID]
                        )  # Change the reference in the dictionary to account for the new ID
                        tmpID = RxnList[RxnList[RxnID].ID].ID
                        del RxnList[RxnID]  # remove the old reaction ID
                        RxnID = tmpID

                    ######### Mass Balance the reaction #########
                    ithRxn = RxnList[RxnID]
                    eq , mb_test = RxnParam2Eq(ithRxn, MetList, MetEquiv)
                    LibIni = WrapRxnSubsProdParam(ithRxn, MetList, MetEquiv)
                    if mb_test != 0:
                        IthRxnMB = mass_balance(eq, RxnID)
                        if IthRxnMB[4]: # If new compounds have to be added to mass balance the reactions, then check if they need to be added to the network as compounds
                            for x in IthRxnMB[4]:
                                if not extra_compound[x[0]] in MetIdent and not extra_compound[x[0]] in MetEquiv:
                                    MetIdent = MetIdent + [extra_compound[x[0]]]
                                    if not x[0] in 'R' and not x[0] in 'X':
                                        extra_url = 'https://www.genome.jp/entry/'+extra_compound[x[0]]
                                        MetList[extra_compound[x[0]]] = compound(extra_url,extra_compound[x[0]],time,EF,specialCompounds)
                                    else:
                                        MetList[extra_compound[x[0]]] = add_extra_compound (x[0],extra_compound, time, EF, specialCompounds)
                    else: # if the reaction cannot be mass balanced all the stoichimetric coef are assumed to be like in the original reaction 
                        IthRxnMB = ([float(x[0]) for x in ithRxn.Substrate()],[float(x[0]) for x in ithRxn.Product()],0,0,0,0,[MetEquiv[x[2]] if x[2] in MetEquiv else x[2] for x in ithRxn.Substrate()], [MetEquiv[x[2]] if x[2] in MetEquiv else x[2] for x in ithRxn.Product()],'','',0)
                    LibEnd = UnwrapRxnSubsProdParam(IthRxnMB, LibIni, IthRxnMB)

                    # Add metabolites and stc coeff to reaction
                    S = list()
                    for x in LibEnd[0]:
                        S.append(
                            [
                                str(LibEnd[0][x][0]),
                                (
                                    "http://www.genome.jp/dbget-bin/www_bget?cpd:"
                                    + LibEnd[0][x][1]
                                ),
                                LibEnd[0][x][1],
                            ]
                        )
                    P = list()
                    for x in LibEnd[1]:
                        P.append(
                            [
                                str(LibEnd[1][x][0]),
                                (
                                    "http://www.genome.jp/dbget-bin/www_bget?cpd:"
                                    + LibEnd[1][x][1]
                                ),
                                LibEnd[1][x][1],
                            ]
                        )
                    S2 = copy.deepcopy(S)
                    P2 = copy.deepcopy(P)
                    RxnList[RxnID].Substrate = lambda: S2
                    RxnList[RxnID].Product = lambda: P2
                    RxnList[RxnID].SetSubstrate = lambda: S2
                    RxnList[RxnID].SetProduct = lambda: P2
                    
                    RxnList[RxnID].subs = S2
                    RxnList[RxnID].prods = P2
                    if not RxnID in PathNameRxn.get(PathName):
                        PathNameRxn[PathName] += RxnID + " "

                    ######### Define New GPR ###########
                    tmpGPR = ()
                    tmpSC = ()
                    for x in RxnList[RxnID].EC():
                        if x not in GPRIdent:
                            GPRIdent = GPRIdent + [x]
                            GPRList[x] = gpr(x, session)
                        try:
                            tmpGPR = tmpGPR + GPRList[x].GprSubcell()[0:2]
                            tmpSC = tmpSC + GPRList[x].GprSubcell()[2:4]
                        except:
                            tmpGPR = tmpGPR + GPRList[x].GprSubcell[0:2]
                            tmpSC = tmpSC + GPRList[x].GprSubcell[2:4]                                
                    # Reorganize S-GPRs and GPRs based on their specific location
                    tmpSC2 = [dict(),dict()]
                    reactio_compartment_list = list(set([x.strip() for x in str([list(x.keys()) for x in tmpSC]).replace('[','').replace(']','').replace('\'','').split(',')]))
                    for x in reactio_compartment_list:
                        xth_tmp_gpr = [y for y in tmpSC if x in y.keys()]
                        tmp_xth_sgpr = ""
                        tmp_xth_gpr = ""
                        for y in range(int(len(xth_tmp_gpr)/2)):
                            if not re.findall('^\[\]$', xth_tmp_gpr[y+y][x]):
                                tmp_xth_sgpr += xth_tmp_gpr[y+y][x]
                            if not re.findall('^\[\]$', xth_tmp_gpr[y+y+1][x]):    
                                tmp_xth_gpr += xth_tmp_gpr[y+y+1][x]
                        tmp_xth_sgpr = str(set(tmp_xth_sgpr.replace('][', '] or [').split(" or "))).replace('\'', '').replace('{', '').replace('}', '').replace(',', ' or')
                        tmp_xth_gpr = str(set(tmp_xth_gpr.replace('][', '] or [').split(" or "))).replace('\'', '').replace('{', '').replace('}', '').replace(',', ' or')
                        tmpSC2[0][x] = tmp_xth_sgpr
                        tmpSC2[1][x] = tmp_xth_gpr

                    RxnList[RxnID].GPR = tmpGPR
                    RxnList[RxnID].Subcel = tmpSC2

                    ######### Expand the annotations based on the cellular location ###########
                    Compartment_CL, rxn_cl, comp_cl = rxnSubcel(
                        RxnList[RxnID],
                        RxnList_CL,
                        MetList_CL,
                        RxnIdent_CL,
                        MetIdent_CL,
                        Compartment_CL,
                        RxnList,
                        MetList,
                        MetEquiv,
                    )
                
                    RxnIdent_CL = RxnIdent_CL + rxn_cl
                    MetIdent_CL = MetIdent_CL + comp_cl

                    
                    print(
                    "Reaction("
                    + str(j + 1)
                    + "/"
                    + str(len(PathList[PathID].Reactions()))
                    + ")_Pathway"
                    + "("
                    + str(i + 1)
                    + "/"
                    + str(len(Path) - 1)
                    + ")"
                    )
            except:
                continue
            j = j + 1
    else:
        print(
            PathName
            + ": not defined in human"
            + "("
            + str(i + 1)
            + "/"
            + str(len(Path) - 1)
            + ")"
        )

    i = i + 1



######### Genes ###########
GeneList = {}
GeneIdent = []
g = 0
while g < len(GPRList):
    if GPRIdent[g] in GPRList.keys() and GPRList[GPRIdent[g]].GprSubcell():
        gene_matches = re.findall(
            "([A-Za-z0-9\-]+)",
            GPRList[GPRIdent[g]].GprSubcell()[1].replace("and", "").replace("or", ""),
        )
        z = 0
        while z < len(gene_matches):
            if not gene_matches[z] in GeneIdent:
                GeneIdent = GeneIdent + [gene_matches[z]]
                GeneList[gene_matches[z]] = gene(gene_matches[z], EnsblDB)
            z = z + 1
    g = g + 1

with open("files/pre_sbml_raw.pk", "wb") as f:
    dill.dump(
        {
            "name": ModName,
            "id": ModID,
            "mets": MetList,
            "mets_cl": MetList_CL,
            "met_equiv": MetEquiv,
            "reactions": RxnList,
            "reactions_cl": RxnList_CL,
            "genes": GeneList,
            "pathways": PathNameRxn,
            "loc": LocVar,
        },
        f,
    )


Compartment_CL = sorted(Compartment_CL)

listOfID = list(CSL_ID.values())  # abbr. id
LipidMasterlistOfID = []

for CSL in Compartment_CL:
    CSL2 = re.sub(r'[^A-Za-z0-9 ]+', '', CSL)
    ModMaster = list(set(LipidMasterlistOfID + listOfID))
    if CSL_ID.get(CSL):
        ID = CSL_ID.get(CSL)
    elif len(re.sub(" $", "", re.sub("^ ", "", CSL)).split(" ")) > 1:
        ID = (CSL2.split(" ")[0][0] + CSL2.split(" ")[1][0]).lower().replace(" ", "")
    else:
        if len(CSL.split(" ")) > 1:
            ID = CSL2[0:3].lower().replace(" ", "")
        else:
            ID = CSL2[0:2].lower().replace(" ", "")
    if not CSL_ID.get(CSL) and ID in ModMaster:
        r = re.compile(ID)
        ID = ID + str(len(list(filter(r.match, ModMaster))) + 1)
    LocVar[CSL] = ""
    LocVar[CSL] += ID
    listOfID.append(ID)
    LipidMasterlistOfID.append(ID)


with open("files/pre_sbml_pos_comp.pk", "wb") as f:
    dill.dump(
        {
            "name": ModName,
            "id": ModID,
            "mets": MetList,
            "mets_cl": MetList_CL,
            "met_equiv": MetEquiv,
            "reactions": RxnList,
            "reactions_cl": RxnList_CL,
            "genes": GeneList,
            "pathways": PathNameRxn,
            "loc": LocVar,
        },
        f,
    )

with open('files/pre_sbml_pos_comp.pk', 'rb') as f:
    data = pickle.load(f)
#with open('files/pre_sbml_pos_comp.pk', 'rb') as f:
#  data = f.read()

model = cobra_reconstruction(
    ModName, ModID, MetList_CL, RxnList_CL, GeneList, PathNameRxn, LocVar, MetEquiv, MetList
)
cobra.io.write_sbml_model(model, Output)



