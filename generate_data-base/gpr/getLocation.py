############## Libraries ##############
import os
import pickle
import random
import re
import time
import urllib
from itertools import zip_longest
from typing import Optional, Union
import urllib.request
import requests
import logging

from functions_generate_database import *


############## Functions ##############
# Option 4: Based on Option 1 using urllib. Less powerful of all but keeps the structure of the output which is very convenient for an easy implementation in the pipeline (https://stackoverflow.com/questions/1096379/how-to-make-urllib2-requests-through-tor-in-python)

LOGGER = logging.getLogger(__name__)
session = setup_biocyc_session()


# class or module with a get() method
HttpGetter = "Getter"


"Load content from html page"


def getHtml(url: str, session: HttpGetter) -> str:
    """
    Loads the content of a web page given an url. it can be done using tor network or not

    Inputs
    ------
    url: url to be loaded in string format
    hide_ip: boolean [0,1] indicating if the tor network has to be activated or not (1 or 0 respectively)

    Output
    ------
    url content as str
    """
    try:
        return HttpGetter.get(url).text
    except Exception as e:
        return ""


def create_dict(gene, e):
    """Create a dictionary of genes."""
    n = ""
    global my_dict
    if "my_dict" not in globals():
        my_dict = {}
    elif gene and e:
        my_dict[gene] = e
    elif gene and not e:
        n = my_dict[gene]
    return my_dict, n


def multiple_replace(dictt, text):
    """Multiple replacement."""
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dictt.keys())))
    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dictt[mo.string[mo.start() : mo.end()]], text)


def getLocation(
    gpr,
    genelist1,
    genelist2,
    impose_locations,
    location_dict_file,
    session: Optional[requests.Session] = None,
):
    """Finds subcellular location.
    Using a SGPR (or GPR), it identifies the corresponding cellular locations and adapts the location-specific SGPRs accordingly.
    i.e.:
        - input GPR: a*1 or b*1 (tipically output from getGPR function)
        - identified cellular location: cytosol, mitochondria
        - identified isoforms in the citosol: a
        - identified isoforms in the mitochondria: b
        - output (cellular location: SGPR): {cytosol: a*1} {mitochondria: b*1}

    Inputs
    ----------
    gpr: url to be loaded in string format
    genelist1: list of gene names (i.e. ['ADH1B', 'ADH6', 'ADH4'])
    genelist2: list of gene biocyc IDs (i.e. ['HS08983', 'HS10600', 'HS06569'])
    impose_locations: boolean. 0: free search. 1: limit the cellular locations to a list of locations and all the locations that do not fit within this list are considered cytosol.
    location_dict_file: file name (".pickle" extension) with the list of compartments

    Output
    -------
    RuleLoc: {cellular location: SGPR}, i.e: {'Cytosol': '[CES2*1]'}
    RuleLoc2: {cellular location: GPR}, i.e: {'Cytosol': '[CES2]'}
    RuleLoc3: {cellular location: GPR with Ensembl IDs}, i.e: {'Cytosol': 'ENSG00000172831'}
    RuleLoc4:  {genes Ensemble ID: genes name}, i.e {'ENSG00000172831': ['CES2']})
    """
    if session is None:
        session = requests
    try:
        gprgpr = ""
        OtherLocations = ["Other locations"]
        gpr2 = gpr
        gpr = re.sub(r"\*[0-9]+", "", gpr)
        urls0 = [
            x.replace("[", "")
            .replace("(", "")
            .replace("]", "")
            .replace(")", "")
            .replace(",", "")
            .replace(" ", "")
            for x in gpr2.split("or")
        ]
        hh = 1
        if impose_locations == 1:
            Var = open(location_dict_file, "rb")  # Endo1b_variables.pkl
            dictt = pickle.load(Var)
            dictt = dict((k.lower(), v.lower()) for k, v in dictt.items())
            hh = 0  # dict and not all
            print(location_dict_file + " is used")
        LocationList = []
        Locations = []
        uniprot = ""
        biocyc = ""
        n = 0
        for n in range(len(urls0)):  # isoforms
            a = urls0[n].split("and")
            m = 0
            d = ""
            for m in range(len(a)):
                p = 1
                b = (
                    "http://www.genome.jp/dbget-bin/www_bget?sp:"
                    + re.sub(r"\*[0-9]+", "", a[m])
                    + "_HUMAN"
                )  # location in genome net human
                # bb = str(getHtml(b, session))  # .decode('utf-8')
                bb = str(urllib.request.urlopen(b).read())
                dd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                if dd:
                    print(dd)
                if not dd:
                    ddd = re.search("(SUBCELLULAR LOCATION:.*)", "")
                    if ddd:
                        dd = [ddd.group(1)]
                if not dd or dd:
                    if genelist1:
                        if isinstance(genelist1, str):
                            genelist11 = re.findall(r"\[(.+?)\]", genelist1)
                            genelist11 = [
                                gene.replace("(", "")
                                .replace(")", "")
                                .replace("[", "")
                                .replace("]", "")
                                for gene in genelist11
                            ]
                        if not isinstance(genelist1, str):
                            genelist11 = genelist1
                        index = [x.upper() for x in genelist11].index(
                            re.sub(r"\*[0-9]+", "", a[m].upper())
                        )
                        b = (
                            "http://biocyc.org/gene?orgid=META&id="
                            + genelist2[index].upper()
                        )  # location in BioCyc human #########################
                        # bb = str(getHtml(b, session))
                        try:
                            bb = str(urllib.request.urlopen(b).read())
                        except urllib.error.HTTPError:
                            bb = ''
                        if not bb:
                            b = (
                                "https://biocyc.org/gene?orgid=HUMAN&id="
                                + genelist2[index].upper()
                            )  # location in BioCyc human #########################
                            #bb = str(getHtml(b, session))
                            try:
                                bb = str(urllib.request.urlopen(b).read())
                            except urllib.error.HTTPError:
                                bb = ''
                            # bb = str(urllib.request.urlopen(b).read())
                        if bb:
                            if re.search("Location", bb, flags=re.DOTALL):
                                ddd1 = re.findall("Locations?.+?Reactions?", bb)[0]
                                ddd = ddd = [
                                    re.sub(r"\\n", "", i)
                                    for i in re.findall(
                                        "([\\\\n+]?[a-z ]+[a-zA-Z() ]+) <a", ddd1
                                    )
                                ]
                                if not ddd:
                                    ddd = [
                                        re.sub(r"\\n", "", i)
                                        for i in re.findall(
                                            "(\\\\[a-zA-Z, ()]+)</", ddd1
                                        )
                                    ]
                                    ddd = [
                                        re.sub("^ ", "", i) for i in ddd[0].split(",")
                                    ]
                                dd.extend(ddd)
                            biocyc = "1"
                            if bb:  # location in Uniprot
                                biocyc = ""
                                if re.search("uniprot/([A-Z0-9]+)", bb):
                                    cc = re.search("uniprot/([A-Z0-9]+)", bb).group(1)
                                    b = (
                                        "https://rest.uniprot.org/uniprotkb/"
                                        + cc
                                        + ".txt"
                                    )
                                    # bb = str(getHtml(b, session))
                                    bb = str(urllib.request.urlopen(b).read())
                                    ddd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                                    ddd = [x for x in ddd if "GO" not in x]
                                    dd.extend(ddd)
                                    uniprot = "1"
                        if not dd:
                            UniProtKB = ""
                            GeneID = re.sub(r"\*[0-9-]+", "", a[m])
                            https = (
                                "https://www.uniprot.org/uniprot/?query="
                                + GeneID
                                + "&sort=score"
                            )
                            # url = str(
                            #     getHtml(https, session)
                            # )  # str(urllib.request.urlopen(https).read())
                            url = str(urllib.request.urlopen(https).read())
                            if re.search(
                                'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                + GeneID
                                + "_HUMAN",
                                url,
                            ):
                                UniProtKB = re.search(
                                    'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                    + GeneID
                                    + "_HUMAN",
                                    url,
                                ).group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE):
                                    UniProtKB = search.group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE)[
                                    0
                                ]:
                                    UniProtKB = search.group(1)
                            cc = UniProtKB
                            b = (
                                "https://www.uniprot.org/uniprot/"
                                + cc
                                + "#subcellular_location"
                            )
                            # bb = str(getHtml(b, session))
                            bb = str(urllib.request.urlopen(b).read())
                            ddd = re.findall(
                                'class="[a-zA-Z_ ]+"><h6>([a-zA-Z ]+)</h6>', bb
                            )
                            if ddd == ["Other locations"]:
                                ddd = re.findall(
                                    'locations*/SL-[0-9]+">([a-zA-Z ]+) </a>', bb
                                )
                            dd.extend(ddd)
                            uniprot = "1"
                dd = [i.strip() for i in dd]
                dd = [d for d in dd if d]
                if dd:
                    ee = [
                        re.sub("}[.,;]", "}@", x, flags=re.DOTALL)
                        .replace("  ", "")
                        .replace("@ ", "")
                        .replace(" {", "{")
                        .replace("SUBCELLULAR LOCATION: ", "")
                        .split("}")[0]
                        for x in dd
                    ]
                    if uniprot or biocyc:
                        ee = re.findall(r"([A-Za-z ]+)", ee[0])
                        uniprot = ""
                    SubUnLoc = []
                    x = 0
                    for x in range(len(dd)):
                        if (
                            dd[x]
                            and not dd[x] in OtherLocations
                            and not "GO" in dd[x]
                            and not "PubMed" in dd[x]
                        ):
                            ee = dd
                            ff = ee[x].split("{")[0]
                            if biocyc:
                                ff = re.findall(r"([A-Za-z0-9 ,;\-\(\)]+)", ff)[0]
                                if re.findall(r",", ff):
                                    ee = ee + ff.split(",")[1:]
                                    ff = ff.split(",")[0]
                                biocyc = ""
                            if re.findall(
                                ":", ff
                            ):  # if there is ":" in the subcellular location, the real location is after :
                                a = ff.split(":")
                                ff = a[1]
                            if hh == 0:
                                ff = ff.lower()
                                ff = multiple_replace(dictt, ff.split("{")[0])
                            if re.findall("dendriti", ff, re.IGNORECASE):
                                ff = "dendrite"
                            if (
                                re.match("membrane", ff, re.IGNORECASE)
                                or re.findall("plasma membrane", ff, re.IGNORECASE)
                                or re.findall("integral component of membrane", ff, re.IGNORECASE)
                            ):
                                ff = "plasma membrane"
                            if re.findall("eroxisom", ff):
                                ff = "peroxisome"
                            if re.findall("itochondri", ff) or re.findall(
                                "mitochondri", ff
                            ):
                                ff = "mitochondria"
                            if re.findall("lysosom", ff) or re.findall("Lysosom", ff):
                                ff = "lysosome"
                            if re.findall("olgi", ff):
                                ff = "golgi apparatus"
                            if re.findall("xtracel", ff):
                                ff = "extracellular"
                            if re.findall("ndoplasm", ff) or re.findall("ndosom", ff):
                                ff = "endoplasmic reticulum"
                            if (
                                re.findall("ytosol", ff)
                                or re.findall("ytoplasm", ff)
                                or re.findall("ntracel", ff)
                            ):
                                ff = "cytosol"
                            if (
                                re.findall("ucle", ff)
                                or re.findall("enter", ff)
                                or re.findall("entro", ff)
                                or re.findall("entri", ff)
                                or re.findall("pindle", ff)
                                or re.findall("rna", ff)
                                or re.findall("dna", ff)
                                or re.findall("smn complex", ff)
                                or re.findall("risc complex", ff)
                                or re.findall("axon", ff)                                
                            ):
                                ff = "nucleus"
                            if hh == 1 and not ff.lower():  # important
                                ff = "cytosol"
                            if hh == 0 and not ff.lower() in list(
                                set(dictt.keys())
                            ):  # important
                                ff = "cytosol"
                            if re.findall("\\\\n", ff) or re.findall(
                                "\.", ff
                            ):  # important
                                ff = "cytosol"
                            if not re.findall("Note=", ff):
                                SubUnLoc = SubUnLoc + [ff]
                                Locations = Locations + [ff]
                    d = sorted(set(SubUnLoc))
                    LocList = [
                        "extracellular",
                        "peroxisome",
                        "mitochondria",
                        "cytosol",
                        "lysosome",
                        "endoplasmic reticulum",
                        "golgi apparatus",
                        "nucleus",
                        "inner mitochondria",
                    ]
                    m = len(
                        a
                    )  # once it is defined a cellular location the process stops because it is assumed that all the subunit of the same complex are in the same place
                if not dd:
                    print(re.sub(r"\*[0-9]+", "", a[m]))
                    m = m + 1
                    if not gpr in gprgpr:
                        gprgpr += gpr + "\n"
                    d = [
                        "cytosol"
                    ]  # by default, if there is not anotated location, the reaction is located into the cytoso
                    p = 0
                    Locations = Locations + [d][0]
            LocationList = LocationList + [d]
            n = n + 1
        Locations = list(set(Locations))
        l = 0
        RuleLoc = {}
        RuleLoc2 = {}
        RuleLoc3 = {}
        RuleLoc4 = {}
        while l < len(Locations):
            L = Locations[l]
            k = 0
            LocGPR = "["
            while k < len(LocationList):
                g = 0
                while g < len(LocationList[k]):
                    if L == LocationList[k][g]:
                        LocGPR = LocGPR + urls0[k] + "] or ["
                    g = g + 1
                k = k + 1
            LocGPR = LocGPR[:-5].replace("and", " and ").replace("  ", " ")
            RuleLoc[Locations[l]] = LocGPR
            RuleLoc2[Locations[l]] = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\[", "", LocGPR2)
            LocGPR2 = re.sub("\]", "", LocGPR2)
            m = LocGPR2.split(" ")[1::2]
            LocGPR2 = LocGPR2.split(" ")[::2]
            GPR6 = create_dict(str(), str())[0]
            for i in range(len(LocGPR2)):
                iiii = ""
                if LocGPR2[i] in GPR6:
                    iiii = create_dict(LocGPR2[i], str())[1]
                else:
                    if LocGPR2[i]:
                        iiii = LocGPR2[i]
                        iiii2 = ""
                        iiii2 = re.search(
                            "gene=([A-Z0-9]+)",
                            str(
                                urllib.request.urlopen(
                                    "https://www.genome.jp/dbget-bin/www_bget?hsa+"
                                    + LocGPR2[i]
                                ).read()
                            ),
                        )  # .group(1)
                        if not iiii2:
                            iiii2 = re.search(
                                "(ENSG[0-9]+)",
                                str(
                                    urllib.request.urlopen(
                                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                                        + LocGPR2[i]
                                    ).read()
                                ),
                            )
                        if iiii2:
                            iiii = iiii2.group(1)
                RuleLoc4[iiii] = list()
                RuleLoc4[iiii].append(LocGPR2[i])
                create_dict(LocGPR2[i], iiii)
                LocGPR2[i] = iiii
            LocGPR2 = " ".join(
                [m + " " + str(n) for m, n in zip_longest(LocGPR2, m, fillvalue="")]
            )[:-1]
            RuleLoc3[Locations[l]] = LocGPR2
            l = l + 1
        if not RuleLoc4 == {"": [""]}:
            print(RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4)
        return (
            RuleLoc,
            RuleLoc2,
            RuleLoc3,
            RuleLoc4,
        )  # , p # p indicates that the location couldn't be determined and citosol has been put instead
    except Exception as e:
        return ""


############## Test Function ##############
# input_pickle_file = "ec_gpr_location.pickle"  # If this file is used (impose_location = 1) all the compartments found that do not fit with the list of compartments defined here will be considered as "cytosol"
# input_variables = open(input_pickle_file, "rb")
# dictt_ini = pickle.load(input_variables)
# dictt_ini_key_list = list(dictt_ini)
# impose_locations = 0
# location_dict_file = "Human_variables.pkl"  # Endo1a_variables.pkl #Human_variables.pkl
#
# t = 0
# for x in dictt_ini_key_list:
#     gpr = dictt_ini[x][3]
#     genelist1 = dictt_ini[x][1]
#     genelist2 = dictt_ini[x][2]
#     get_loc = getLocation(
#         gpr, genelist1, genelist2, impose_locations, location_dict_file
#     )
#     print(t, get_loc)
#     t = t + 1
#
#
# ############## Function development ##############
# os.chdir(
#     "/home/igor/Documentos/Full-Human-GEM/ENDO-MODS_SGPR_Improvement_2022-10-21/DEV_Folder"
# )
#
# # Extract input data from pickle file
# input_pickle_file = "ec_gpr_location.pickle"
# input_variables = open(input_pickle_file, "rb")
# dictt_ini = pickle.load(input_variables)
# dictt_ini_key_list = list(dictt_ini)
#
# # Explore what is the best example to develope the function
# t = 0
# for x in dictt_ini_key_list:
#     print(t, dictt_ini[x][0], dictt_ini[x][1], dictt_ini[x][2])
#     t = t + 1
#     0
#
# v = 778
# gpr = dictt_ini[dictt_ini_key_list[v]][3]
# genelist1 = dictt_ini[dictt_ini_key_list[v]][1]
# genelist2 = dictt_ini[dictt_ini_key_list[v]][2]
# impose_locations = 0
# location_dict_file = "Human_variables.pkl"  # Endo1a_variables.pkl #Human_variables.pkl
#
# getLocation(gpr, genelist1, genelist2, impose_locations, location_dict_file)

# Final:
""""Finds subcellular location"""


def getLocation_final(gpr, genelist1, genelist2, impose_locations, location_dict_file):
    try:
        hide_ip = 0
        gprgpr = ""
        OtherLocations = ["Other locations"]
        gpr2 = gpr
        gpr = re.sub(r"\*[0-9]+", "", gpr)
        urls0 = [
            x.replace("[", "")
            .replace("(", "")
            .replace("]", "")
            .replace(")", "")
            .replace(",", "")
            .replace(" ", "")
            for x in gpr2.split("or")
        ]
        hh = 1
        if impose_locations == 1:
            Var = open(location_dict_file, "rb")  # Endo1b_variables.pkl
            dictt = pickle.load(Var)
            hh = 0  # dict and not all
            print(location_dict_file + " is used")
        LocationList = []
        Locations = []
        uniprot = ""
        biocyc = ""
        n = 0
        for n in range(len(urls0)):  # isoforms
            a = urls0[n].split("and")
            m = 0
            d = ""
            for m in range(len(a)):
                p = 1
                b = (
                    "http://www.genome.jp/dbget-bin/www_bget?sp:"
                    + re.sub(r"\*[0-9]+", "", a[m])
                    + "_HUMAN"
                )  # location in genome net human
                print("A")
                bb = str(getHtml(b, session))  # .decode('utf-8')
                dd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                if dd:
                    print(dd)
                if not dd:
                    ddd = re.search("(SUBCELLULAR LOCATION:.*)", "")
                    if ddd:
                        dd = [ddd.group(1)]
                if not dd or dd:
                    if genelist1:
                        if isinstance(genelist1, str):
                            genelist11 = re.findall("\[(.+?)\]", genelist1)
                            genelist11 = [
                                gene.replace("(", "")
                                .replace(")", "")
                                .replace("[", "")
                                .replace("]", "")
                                for gene in genelist11
                            ]
                        if not isinstance(genelist1, str):
                            genelist11 = genelist1
                        index = [x.upper() for x in genelist11].index(
                            re.sub(r"\*[0-9]+", "", a[m].upper())
                        )
                        b = (
                            "http://biocyc.org/gene?orgid=META&id="
                            + genelist2[index].upper()
                        )  # location in BioCyc human #########################
                        print("B")
                        bb = str(getHtml(b, session))
                        if not bb:
                            b = (
                                "https://biocyc.org/gene?orgid=HUMAN&id="
                                + genelist2[index].upper()
                            )  # location in BioCyc human #########################
                            print("C")
                            bb = str(getHtml(b, session))
                        if bb:
                            if re.search("Location", bb, flags=re.DOTALL):
                                ddd1 = re.findall("Locations?.+?Reactions?", bb)[0]
                                ddd = ddd = [
                                    re.sub(r"\\n", "", i)
                                    for i in re.findall(
                                        "([\\\\n+]?[a-z ]+[a-zA-Z() ]+) <a", ddd1
                                    )
                                ]
                                if not ddd:
                                    ddd = [
                                        re.sub(r"\\n", "", i)
                                        for i in re.findall(
                                            "(\\\\[a-zA-Z, ()]+)</", ddd1
                                        )
                                    ]
                                    ddd = [
                                        re.sub("^ ", "", i) for i in ddd[0].split(",")
                                    ]
                                dd.extend(ddd)
                            biocyc = "1"
                            if bb:  # location in Uniprot
                                biocyc = ""
                                if re.search("uniprot/([A-Z0-9]+)", bb):
                                    cc = re.search("uniprot/([A-Z0-9]+)", bb).group(1)
                                    b = (
                                        "https://rest.uniprot.org/uniprotkb/"
                                        + cc
                                        + ".txt"
                                    )
                                    print("D")
                                    bb = str(getHtml(b, session))
                                    ddd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                                    ddd = [x for x in ddd if not "GO" in x]
                                    dd.extend(ddd)
                                    uniprot = "1"
                        if not dd:
                            UniProtKB = ""
                            GeneID = re.sub(r"\*[0-9-]+", "", a[m])
                            https = (
                                "https://www.uniprot.org/uniprot/?query="
                                + GeneID
                                + "&sort=score"
                            )
                            print("E")
                            url = str(
                                getHtml(https, session)
                            )  # str(urllib.request.urlopen(https).read())
                            if re.search(
                                'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                + GeneID
                                + "_HUMAN",
                                url,
                            ):
                                UniProtKB = re.search(
                                    'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                    + GeneID
                                    + "_HUMAN",
                                    url,
                                ).group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE):
                                    UniProtKB = search.group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE)[
                                    0
                                ]:
                                    UniProtKB = search.group(1)
                            cc = UniProtKB
                            b = (
                                "https://www.uniprot.org/uniprot/"
                                + cc
                                + "#subcellular_location"
                            )
                            print("F")
                            bb = str(getHtml(b, session))
                            ddd = re.findall(
                                'class="[a-zA-Z_ ]+"><h6>([a-zA-Z ]+)</h6>', bb
                            )
                            if ddd == ["Other locations"]:
                                ddd = re.findall(
                                    'locations*/SL-[0-9]+">([a-zA-Z ]+) </a>', bb
                                )
                            dd.extend(ddd)
                            uniprot = "1"
                dd = [i.strip() for i in dd]
                dd = [d for d in dd if d]
                if dd:
                    ee = [
                        re.sub("}[.,;]", "}@", x, flags=re.DOTALL)
                        .replace("  ", "")
                        .replace("@ ", "")
                        .replace(" {", "{")
                        .replace("SUBCELLULAR LOCATION: ", "")
                        .split("}")[0]
                        for x in dd
                    ]
                    if uniprot or biocyc:
                        ee = re.findall(r"([A-Za-z ]+)", ee[0])
                        uniprot = ""
                    SubUnLoc = []
                    x = 0
                    for x in range(len(dd)):
                        if (
                            dd[x]
                            and not dd[x] in OtherLocations
                            and not "GO" in dd[x]
                            and not "PubMed" in dd[x]
                        ):
                            ee = dd
                            ff = ee[x].split("{")[0]
                            if biocyc:
                                ff = re.findall(r"([A-Za-z0-9 ,;\-\(\)]+)", ff)[0]
                                if re.findall(r",", ff):
                                    ee = ee + ff.split(",")[1:]
                                    ff = ff.split(",")[0]
                                biocyc = ""
                            if re.findall(
                                ":", ff
                            ):  # if there is ":" in the subcellular location, the real location is after :
                                a = ff.split(":")
                                ff = a[1]
                            if hh == 0:
                                ff = ff.lower()
                                ff = multiple_replace(dictt, ff.split("{")[0])
                            if re.findall("Dendriti", ff, re.IGNORECASE):
                                ff = "Dendrite"
                            if (
                                re.match("membrane", ff, re.IGNORECASE)
                                or re.findall("plasma membrane", ff, re.IGNORECASE)
                                or re.findall(
                                    "integral component of membrane", ff, re.IGNORECASE
                                )
                            ):
                                ff = "Plasma membrane"
                            if re.findall("eroxisom", ff):
                                ff = "Peroxisome"
                            if re.findall("itochondri", ff) or re.findall(
                                "mitochondri", ff
                            ):
                                ff = "Mitochondria"
                            if re.findall("lysosom", ff) or re.findall("Lysosom", ff):
                                ff = "Lysosome"
                            if re.findall("olgi", ff):
                                ff = "Golgi apparatus"
                            if re.findall("xtracel", ff):
                                ff = "Extracellular"
                            if re.findall("ndoplasm", ff) or re.findall("ndosom", ff):
                                ff = "Endoplasmic reticulum"
                            if (
                                re.findall("ytosol", ff)
                                or re.findall("ytoplasm", ff)
                                or re.findall("ntracel", ff)
                            ):
                                ff = "cytosol"
                            if (
                                re.findall("ucle", ff)
                                or re.findall("enter", ff)
                                or re.findall("entro", ff)
                                or re.findall("entri", ff)
                                or re.findall("pindle", ff)
                                or re.findall("RNA", ff)
                                or re.findall("DNA", ff)
                                or re.findall("SMN complex", ff)
                                or re.findall("RISC complex", ff)
                                or re.findall("axon", ff)
                                or re.findall("Axon", ff)
                            ):
                                ff = "Nucleus"
                            if hh == 1 and not ff.lower():  # important
                                ff = "cytosol"
                            if hh == 0 and not ff.lower() in list(
                                set(dictt.keys())
                            ):  # important
                                ff = "cytosol"
                            if re.findall("\\\\n", ff) or re.findall(
                                "\.", ff
                            ):  # important
                                ff = "cytosol"
                            if not re.findall("Note=", ff):
                                SubUnLoc = SubUnLoc + [ff]
                                Locations = Locations + [ff]
                    d = sorted(set(SubUnLoc))
                    LocList = [
                        "extracellular",
                        "peroxisome",
                        "mitochondria",
                        "cytosol",
                        "lysosome",
                        "endoplasmic reticulum",
                        "golgi apparatus",
                        "nucleus",
                        "inner mitochondria",
                    ]
                    m = len(
                        a
                    )  # once it is defined a cellular location the process stops because is assumed that all the subunit of the same complex are in the same place
                if not dd:
                    print(re.sub(r"\*[0-9]+", "", a[m]))
                    m = m + 1
                    if not gpr in gprgpr:
                        gprgpr += gpr + "\n"
                    d = [
                        "cytosol"
                    ]  # by default, if there is not anotated location, the reaction is located into the cytoso
                    p = 0
                    Locations = Locations + [d][0]
            LocationList = LocationList + [d]
            n = n + 1
        Locations = list(set(Locations))
        l = 0
        RuleLoc = {}
        RuleLoc2 = {}
        RuleLoc3 = {}
        RuleLoc4 = {}
        while l < len(Locations):
            L = Locations[l]
            k = 0
            LocGPR = "["
            while k < len(LocationList):
                g = 0
                while g < len(LocationList[k]):
                    if L == LocationList[k][g]:
                        LocGPR = LocGPR + urls0[k] + "] or ["
                    g = g + 1
                k = k + 1
            LocGPR = LocGPR[:-5].replace("and", " and ").replace("  ", " ")
            RuleLoc[Locations[l]] = LocGPR
            RuleLoc2[Locations[l]] = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\[", "", LocGPR2)
            LocGPR2 = re.sub("\]", "", LocGPR2)
            m = LocGPR2.split(" ")[1::2]
            LocGPR2 = LocGPR2.split(" ")[::2]
            GPR6 = create_dict(str(), str())[0]
            for i in range(len(LocGPR2)):
                iiii = ""
                if LocGPR2[i] in GPR6:
                    iiii = create_dict(LocGPR2[i], str())[1]
                else:
                    if LocGPR2[i]:
                        iiii = LocGPR2[i]
                        iiii2 = ""
                        iiii2 = re.search(
                            "gene=([A-Z0-9]+)",
                            str(
                                urllib.request.urlopen(
                                    "https://www.genome.jp/dbget-bin/www_bget?hsa+"
                                    + LocGPR2[i]
                                ).read()
                            ),
                        )  # .group(1)
                        if not iiii2:
                            iiii2 = re.search(
                                "(ENSG[0-9]+)",
                                str(
                                    urllib.request.urlopen(
                                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                                        + LocGPR2[i]
                                    ).read()
                                ),
                            )
                        if iiii2:
                            iiii = iiii2.group(1)
                RuleLoc4[iiii] = list()
                RuleLoc4[iiii].append(LocGPR2[i])
                create_dict(LocGPR2[i], iiii)
                LocGPR2[i] = iiii
            LocGPR2 = " ".join(
                [m + " " + str(n) for m, n in zip_longest(LocGPR2, m, fillvalue="")]
            )[:-1]
            RuleLoc3[Locations[l]] = LocGPR2
            l = l + 1
        if not RuleLoc4 == {"": [""]}:
            print(RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4)
        return (
            RuleLoc,
            RuleLoc2,
            RuleLoc3,
            RuleLoc4,
        )  # , p # p indicates that the location couldn't be determined and citosol has been put instead
    except Exception as e:
        return ""


# Original:
def getLocation_original(gpr, genelist1, genelist2, time):
    #    if genelist2: print(genelist2)
    #    print(gpr,genelist1,genelist2,time)
    # print(genelist2)
    gprgpr = ""
    ppList = list()
    OtherLocations = ["Other locations"]
    try:
        # print(genelist1)
        gpr2 = gpr
        gpr = re.sub(r"\*[0-9]+", "", gpr)
        urls0 = [
            x.replace("[", "")
            .replace("(", "")
            .replace("]", "")
            .replace(")", "")
            .replace(",", "")
            .replace(" ", "")
            for x in gpr2.split("or")
        ]
        # print(0, urls0)
        #        dict = {"Cell" : "Cytosol", "Membrane" : "Cytosol", "Lipid-anchor" : "Cytosol", "Cytosol membrane" : "Cytosol", "Multi-pass membrane protein" : "Cytosol", "Single-pass membrane protein" : "Cytosol","terminal bouton":"Cytosol","catalytic complex":"Cytosol", "endocytic vesicle":"Cytosol", "presynapse":"Cytosol","neuron projection":"Cytosol","synapse":"Cytosol","glutamatergic synapse":"Cytosol","lamellipodium":"Cytosol","postsynapse":"Cytosol","transport vesicle":"Cytosol","AMPA glutamate receptor complex":"Cytosol","caveola":"Cytosol","excitatory synapse":"Cytosol", " Mitochondrion inner membrane" : "Mitochondrion", "Intermembrane side" : "Mitochondrion", "Endoplasmic reticulum membrane" : "Endoplasmic"} # expand this dict
        Path = "pkl/Human_variables.pkl"  # Endo1a_variables.pkl #Human_variables.pkl
        Var = open(Path, "rb")  # Endo1b_variables.pkl
        hh = 1  # dict and not all
        dictt = pickle.load(Var)
        if hh == 0:
            print(Path + " is used")
            # print(dictt.keys())
        LocationList = []
        Locations = []
        uniprot = ""
        biocyc = ""
        n = 0
        # if urls0: print(urls0)
        for n in range(len(urls0)):  # isoforms
            # print(n)
            # print(urls0)
            # print(urls0[n])
            a = urls0[n].split("and")
            m = 0
            d = ""
            # while m < len(a): # subunits
            for m in range(len(a)):
                # print str(n)+'_'+str(m)+'_00'
                p = 1
                b = (
                    "http://www.genome.jp/dbget-bin/www_bget?sp:"
                    + re.sub(r"\*[0-9]+", "", a[m])
                    + "_HUMAN"
                )  # location in genome net human
                # print(b)
                bb = str(getHtml(b, session))  # .decode('utf-8')
                # cc = bb.replace("\nCC","").replace("-!-","\n")
                # dd = re.findall(r"GO:[0-9]+.*;.*C:(.*);",cc)
                dd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                # print str(n)+'_'+str(m)+'_A'
                # print(dd,1000)
                if dd:
                    print(b)
                    print(dd)
                if not dd:
                    ddd = re.search("(SUBCELLULAR LOCATION:.*)", "")
                    if ddd:
                        dd = [ddd.group(1)]
                        # print(b)
                        # print(dd)
                        # print str(n)+'_'+str(m)+'_A2'
                if not dd or dd:
                    # print(77777)
                    # b = 'http://www.genome.jp/dbget-bin/www_bget?sp:'+re.sub(r"\*[0-9]+" , "", a[m])+'_MOUSE' #location in genome net mouse
                    # print(b)
                    # bb =  str(getHtml(b,time))
                    # cc = bb.replace("\nCC","").replace("-!-","\n")
                    # dd = re.findall(r"GO:[0-9]+.*;.*C:(.*);",cc) # GANC
                    # if dd:
                    #    print(b)
                    #    print('mus')
                    # print str(n)+'_'+str(m)+'_B'
                    # if not dd:
                    #    ddd = re.search('(SUBCELLULAR LOCATION:.*)',cc)
                    #    if ddd:
                    #            dd = [ddd.group(1)]
                    #        ddd = re.search('(SUBCELLULAR LOCATION:.*)',cc)
                    #        print(b)
                    # print(dd)
                    # print(ddd)
                    # print str(n)+'_'+str(m)+'_B2'
                    # dd =[]
                    if genelist1:
                        # print(a[m].upper())
                        # if not type(genelist1) == '''<class 'str'>''': print(type(genelist1))
                        if isinstance(genelist1, str):
                            genelist11 = re.findall("\[(.+?)\]", genelist1)
                            genelist11 = [
                                gene.replace("(", "")
                                .replace(")", "")
                                .replace("[", "")
                                .replace("]", "")
                                for gene in genelist11
                            ]
                        if not isinstance(genelist1, str):
                            genelist11 = genelist1
                        index = [x.upper() for x in genelist11].index(
                            re.sub(r"\*[0-9]+", "", a[m].upper())
                        )
                        # print(index)
                        # print(genelist2)

                        # print(genelist11[index].upper())
                        b = (
                            "http://biocyc.org/gene?orgid=META&id="
                            + genelist2[index].upper()
                        )  # location in BioCyc human #########################
                        # print(b)
                        # print(genelist11[index].upper(), re.sub(r"\*[0-9]+" , "", a[m].upper()))
                        # if genelist11[index].upper() != re.sub(r"\*[0-9]+" , "", a[m].upper()):    # !!!????
                        # if b:
                        bb = str(getHtml(b, session))
                        if bb:
                            # print(88)
                            # print(d)
                            if re.search("Location", bb, flags=re.DOTALL):
                                print(b)
                                ddd1 = re.findall("Locations?.+?Reactions?", bb)[0]
                                ddd = ddd = [
                                    re.sub(r"\\n", "", i)
                                    for i in re.findall(
                                        "([\\\\n+]?[a-z ]+[a-zA-Z() ]+) <a", ddd1
                                    )
                                ]
                                if not ddd:
                                    ddd = [
                                        re.sub(r"\\n", "", i)
                                        for i in re.findall(
                                            "(\\\\[a-zA-Z, ()]+)</", ddd1
                                        )
                                    ]
                                    # print(ddd)
                                    ddd = [
                                        re.sub("^ ", "", i) for i in ddd[0].split(",")
                                    ]
                                print(ddd)
                                dd.extend(ddd)
                            #                            bb = bytes(bb)
                            #                            cc = cc.encode('utf-8')
                            # cc = str(bb).replace(r"\n",r"").replace(r"<td align=RIGHT valign=TOP class=\"label\">",r"\n<td align=RIGHT valign=TOP class=\"label\"")
                            #                            dd = re.search('<td align=RIGHT valign=TOP class=\"label\"[\S\s]+Location(.*)',cc)
                            # dd = re.findall(r'<td align=RIGHT valign=TOP class=\"label\"[\S\s]+Location(.*)',cc)
                            # if dd:
                            # print(b)
                            #    #print(dd)
                            # dd=''
                            biocyc = "1"
                            # print str(n)+'_'+str(m)+'_D'
                            if bb:  # location in Uniprot
                                biocyc = ""
                                # cc = str(bb).replace(r"\n",r"").replace(r"http://www.uniprot.org",r"\nhttp://www.uniprot.org").replace(r"nbsp",r"nbsp\n")
                                # cc2 = re.findall(r'(http://www.uniprot.org.*)\">',cc)
                                # print(cc2)
                                # if cc2:
                                #    cc3 =  getHtml(cc2[0],time)
                                #    cc4 = cc3.replace("\n","").replace("</span>Subcellular location","\n</span>Subcellular location").replace("&#xd;","\n")
                                #                                #    dd = re.search('</span>Subcellular location(.*)',cc4)
                                #    dd = re.findall('</span>Subcellular location(.*)',cc4)
                                #        if dd:
                                #            print(cc2[0])
                                #            #print(dd)
                                if re.search("uniprot/([A-Z0-9]+)", bb):
                                    cc = re.search("uniprot/([A-Z0-9]+)", bb).group(1)
                                    b = (
                                        "https://rest.uniprot.org/uniprotkb/"
                                        + cc
                                        + ".txt"
                                    )
                                    # print(b)
                                    bb = str(getHtml(b, session))
                                    ddd = re.findall("GO:[0-9]+.+?C:(.+?);", bb)
                                    # print(ddd)
                                    ddd = [x for x in ddd if not "GO" in x]
                                    dd.extend(ddd)
                                    if ddd:
                                        print(b)
                                        print(ddd)
                                    uniprot = "1"
                                    # print str(n)+'_'+str(m)+'_E'
                        if not dd:
                            UniProtKB = ""
                            GeneID = re.sub(r"\*[0-9-]+", "", a[m])
                            # print(GeneID)
                            # print(GeneID)
                            https = (
                                "https://www.uniprot.org/uniprot/?query="
                                + GeneID
                                + "&sort=score"
                            )
                            url = str(urllib.request.urlopen(https).read())
                            # print(https)
                            if re.search(
                                'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                + GeneID
                                + "_HUMAN",
                                url,
                            ):
                                # print(0)
                                UniProtKB = re.search(
                                    'uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>'
                                    + GeneID
                                    + "_HUMAN",
                                    url,
                                ).group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                # print(1)
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_HUMAN.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                # print(search.group())
                                # print(search.group(2))
                                if re.findall(GeneID, search.group(2), re.IGNORECASE):
                                    # print(2)
                                    UniProtKB = search.group(1)
                            elif re.search(
                                '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9-]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                url,
                                re.IGNORECASE,
                            ):
                                # print(3)
                                search = re.search(
                                    '<a href="\/uniprot\/([A-Z0-9-]+)">[A-Z0-9]+<\/a><\/td><td>[A-Z0-9-]+_MOUSE.+?<div class="gene-names"><span class="shortName">(<strong>[a-zA-Z0-9-]+<\/strong>[a-zA-Z0-9-, ]*)<\/span><\/div><\/td><td>?',
                                    url,
                                    re.IGNORECASE,
                                )
                                if re.findall(GeneID, search.group(2), re.IGNORECASE)[
                                    0
                                ]:
                                    # print(4)
                                    UniProtKB = search.group(1)
                            cc = UniProtKB
                            # print(cc)
                            b = (
                                "https://www.uniprot.org/uniprot/"
                                + cc
                                + "#subcellular_location"
                            )
                            # print(b)
                            bb = str(getHtml(b, session))
                            ddd = re.findall(
                                'class="[a-zA-Z_ ]+"><h6>([a-zA-Z ]+)</h6>', bb
                            )
                            # print(ddd)
                            # if ddd:
                            # print(b)
                            if ddd == ["Other locations"]:
                                ddd = re.findall(
                                    'locations*/SL-[0-9]+">([a-zA-Z ]+) </a>', bb
                                )
                                # print(ddd)
                            if GeneID == "Uox":  # uniprot
                                ddd = ["Peroxisome", "Mitochondrion"]  # mouse
                            if GeneID == "NME1-NME2":  # uniprot
                                ddd = ["cytosol"]
                            if GeneID == "CKMT1A" or GeneID == "CKMT1B":  # uniprot
                                ddd = ["mitochondrion"]
                            if GeneID == "TRMT11":  # uniprot
                                ddd = ["cytosol"]
                            dd.extend(ddd)
                            uniprot = "1"
                dd = [i.strip() for i in dd]
                dd = [d for d in dd if d]
                # print(dd)
                if dd:
                    # print(dd)
                    # for i in dd:
                    #    print(i)
                    # print('he')
                    ee = [
                        re.sub("}[.,;]", "}@", x, flags=re.DOTALL)
                        .replace("  ", "")
                        .replace("@ ", "")
                        .replace(" {", "{")
                        .replace("SUBCELLULAR LOCATION: ", "")
                        .split("}")[0]
                        for x in dd
                    ]
                    if uniprot or biocyc:
                        ee = re.findall(r"([A-Za-z ]+)", ee[0])  # ?????!!!!!
                        uniprot = ""
                    SubUnLoc = []
                    x = 0
                    for x in range(len(dd)):
                        if (
                            dd[x]
                            and not dd[x] in OtherLocations
                            and not "GO" in dd[x]
                            and not "PubMed" in dd[x]
                        ):  # ['3.5.1.6']
                            ee = dd
                            ff = ee[x].split("{")[0]
                            if biocyc:
                                ff = re.findall(r"([A-Za-z0-9 ,;\-\(\)]+)", ff)[0]
                                # print(ff)
                                if re.findall(r",", ff):
                                    ee = ee + ff.split(",")[1:]
                                    ff = ff.split(",")[0]
                                biocyc = ""
                            if re.findall(
                                ":", ff
                            ):  # if there is ":" in the subcellular location, the real location is after :
                                a = ff.split(":")
                                ff = a[1]
                            if hh == 0:
                                ff = ff.lower()
                                # print(ff)
                                ff = multiple_replace(
                                    dictt, ff.split("{")[0]
                                )  ### !!!!!!??????
                            # print(ff)
                            # print(dictt)
                            # print(ff, dictt.get(ff))
                            if re.findall("Dendriti", ff, re.IGNORECASE):
                                ff = "Dendrite"
                            if (
                                re.match("membrane", ff, re.IGNORECASE)
                                or re.findall("plasma membrane", ff, re.IGNORECASE)
                                or re.findall(
                                    "integral component of membrane", ff, re.IGNORECASE
                                )
                            ):
                                ff = "Plasma membrane"
                            if re.findall("eroxisom", ff):
                                ff = "Peroxisome"
                            if re.findall("itochondri", ff) or re.findall(
                                "mitochondri", ff
                            ):
                                ff = "Mitochondria"
                            if re.findall("lysosom", ff) or re.findall("Lysosom", ff):
                                ff = "Lysosome"
                            if re.findall("olgi", ff):
                                ff = "Golgi apparatus"
                            if re.findall("xtracel", ff):
                                ff = "Extracellular"
                            if re.findall("ndoplasm", ff) or re.findall("ndosom", ff):
                                ff = "Endoplasmic reticulum"
                            if (
                                re.findall("ytosol", ff)
                                or re.findall("ytoplasm", ff)
                                or re.findall("ntracel", ff)
                            ):
                                ff = "cytosol"
                            if (
                                re.findall("ucle", ff)
                                or re.findall("enter", ff)
                                or re.findall("entro", ff)
                                or re.findall("entri", ff)
                                or re.findall("pindle", ff)
                                or re.findall("RNA", ff)
                                or re.findall("DNA", ff)
                                or re.findall("SMN complex", ff)
                                or re.findall("RISC complex", ff)
                                or re.findall("axon", ff)
                                or re.findall("Axon", ff)
                            ):
                                ff = "Nucleus"
                            # print(ff,list(set(dictt.values())))
                            if hh == 1 and not ff.lower() in list(
                                set(dictt.keys())
                            ):  # important
                                print(1236, ff)
                                # print(sorted(set(dd)))
                                ff = "cytosol"
                            # if hh == 0 and not ff in list(set(dictt.values())): #do not know if this is nessecary
                            #    print(1240, ff)
                            #    #print(sorted(set(dd)))
                            #    ff = 'Cytosol'
                            if re.findall("\\\\n", ff) or re.findall(
                                "\.", ff
                            ):  # important
                                # print(ff)
                                # print(sorted(set(dd)))
                                ff = "Cytosol"
                                print(ff)
                                print()
                            # if hh == 1 and not ff in list(set(dictt.keys())):
                            #    print(ff)
                            #    print(sorted(set(dd)))
                            if not re.findall("Note=", ff):
                                SubUnLoc = SubUnLoc + [ff]
                                Locations = Locations + [ff]

                        # x += 1
                    d = sorted(set(SubUnLoc))
                    # print(d)
                    LocList = [
                        "Extracellular",
                        "Peroxisome",
                        "Mitochondria",
                        "cytosol",
                        "Lysosome",
                        "Endoplasmic reticulum",
                        "Golgi apparatus",
                        "Nucleus",
                        "Inner mitochondria",
                    ]
                    # print([x for x in d if x not in LocList])
                    m = len(
                        a
                    )  # once it is defined a cellular location the process stops because is assumed that all the subunit of the same complex are in the same place
                if not dd:
                    print(re.sub(r"\*[0-9]+", "", a[m]))
                    m = m + 1
                    # if not d and gpr or genelist1 or genelist2:
                    if not gpr in gprgpr:
                        # print()
                        # print('No location was found')
                        # print(gpr, genelist1, genelist2)
                        # print(genelist1,)
                        # print(genelist2)
                        # print()
                        gprgpr += gpr + "\n"
                    d = [
                        "cytosol"
                    ]  # by default, if there is not anotated location, the reaction is located into the cytoso
                    # d = ['']
                    p = 0
                    # ppList.extend(p)
                    Locations = Locations + [d][0]
                # print 'not dd'
            LocationList = LocationList + [d]
            # print(LocationList)
            # print str(n)+'_'+str(m)+'_000'
            n = n + 1
        Locations = list(set(Locations))
        l = 0
        RuleLoc = {}
        RuleLoc2 = {}
        RuleLoc3 = {}
        RuleLoc4 = {}
        # print()
        while l < len(Locations):  # and p != 1:
            L = Locations[l]
            k = 0
            LocGPR = "["
            while k < len(LocationList):
                g = 0
                while g < len(LocationList[k]):
                    # print(0, LocationList)
                    # print(1, urls0)
                    if L == LocationList[k][g]:
                        LocGPR = LocGPR + urls0[k] + "] or ["
                    g = g + 1
                k = k + 1
            LocGPR = LocGPR[:-5].replace("and", " and ").replace("  ", " ")
            RuleLoc[Locations[l]] = LocGPR
            RuleLoc2[Locations[l]] = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\*[0-9]+", "", LocGPR)
            LocGPR2 = re.sub("\[", "", LocGPR2)
            LocGPR2 = re.sub("\]", "", LocGPR2)
            # print(LocGPR2)
            m = LocGPR2.split(" ")[1::2]
            LocGPR2 = LocGPR2.split(" ")[::2]
            GPR6 = create_dict(str(), str())[0]
            # print(0, LocGPR2)
            # print(0, LocGPR2)
            for i in range(len(LocGPR2)):
                # print(0, LocGPR2[i])
                # RuleLoc4[LocGPR2[i]] = list()
                # print(RuleLoc4)
                # print(1, LocGPR2[i])
                # print(LocGPR2[i], 0)
                iiii = ""
                if LocGPR2[i] in GPR6:
                    iiii = create_dict(LocGPR2[i], str())[1]
                    # print(20, iiii)
                else:
                    if LocGPR2[i]:
                        # print('æææ',LocGPR2[i])
                        iiii = LocGPR2[i]
                        iiii2 = ""
                        # print(LocGPR2[i], 1)
                        iiii2 = re.search(
                            "gene=([A-Z0-9]+)",
                            str(
                                urllib.request.urlopen(
                                    "https://www.genome.jp/dbget-bin/www_bget?hsa+"
                                    + LocGPR2[i]
                                ).read()
                            ),
                        )  # .group(1)
                        if not iiii2:
                            iiii2 = re.search(
                                "(ENSG[0-9]+)",
                                str(
                                    urllib.request.urlopen(
                                        "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
                                        + LocGPR2[i]
                                    ).read()
                                ),
                            )
                        if iiii2:
                            iiii = iiii2.group(1)
                            # print(iiii)
                # print(30,  RuleLoc4)
                RuleLoc4[iiii] = list()
                RuleLoc4[iiii].append(LocGPR2[i])
                create_dict(LocGPR2[i], iiii)
                LocGPR2[i] = iiii
                # print(0, LocGPR2[i])
                # print(LocGPR2)
            # print(10, LocGPR2)
            LocGPR2 = " ".join(
                [m + " " + str(n) for m, n in zip_longest(LocGPR2, m, fillvalue="")]
            )[:-1]
            # print(0, LocGPR2)
            RuleLoc3[Locations[l]] = LocGPR2
            # print(3, RuleLoc3)
            l = l + 1
        # print 'FIN'
        # if not RuleLoc3:RuleLoc3= {'1': '1'}
        # print()
        # print(type(RuleLoc3))
        if not RuleLoc4 == {"": [""]}:
            print()
            print(RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4)
            print()
        return (
            RuleLoc,
            RuleLoc2,
            RuleLoc3,
            RuleLoc4,
            p,
        )  # p only with programa_3_2  just a single number and not a list because of line 1256 ca                     m = len(a) # once it is defined a cellular location the process stops because is assumed that all the subunit of the same complex are in the same place
        # print(0, RuleLoc , RuleLoc2)
    except Exception as e:
        # print(traceback.format_exc())
        # RuleLoc = {}
        # RuleLoc2 = {}
        # RuleLoc['Cytosol'] = ''
        # RuleLoc2['Cytosol'] = ''
        # print('Exception "'+str(e)+'" in getLocation with following arguments:')
        # print('gpr: '+str(gpr))
        # print('genelist1: '+str(genelist1))
        # print('genelist2: '+str(genelist2))
        # print(RuleLoc , RuleLoc2)
        return RuleLoc, RuleLoc2, RuleLoc3, RuleLoc4
