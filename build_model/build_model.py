# -*- coding: utf-8 -*-

from cobra.io import read_sbml_model,write_sbml_model
import cobra
import os
import urllib.request, urllib.error, urllib.parse
import re
import urllib.request, urllib.parse, urllib.error
import requests
import copy
import time
import traceback
import itertools
import pubchempy as pcp
import string
import pickle 
from collections import defaultdict
from itertools import zip_longest
from function_build_model import * 
from equations_build_model import *
from cobra import Model, Reaction, Metabolite
from collections import ChainMap
import dill



from gpr.auth_gpr import getGPR, setup_biocyc_session
# 👇 import the addtional function
from gpr.getLocation import getLocation


session = setup_biocyc_session()
location_pkl_file = 'files/bb.pickle'#Human-variables.pkl'



model = read_sbml_model('models/THG-beta1.1.xml')



variablesFile = 'files/variablesFile.pkl' 
Output = 'model_bb.xml'  
Output2 = 'Rxn2Fix'
Output3 = 'THG-beta2-tmp.xml'
Output4 = 'THG-beta2.xml' 
Output5 = 'THG-beta1.xml' 

    
handle = open(Output2, 'w')


c = compartment_file_to_dict() 

variables = defaultdict(list)


listA,listB,ComptoID2 = list(c.values()), [], {}

n = 1 
eList = []
RxnList = []

RxnPath = defaultdict(list)
for p in model.groups:
    for r in p.members:
        RxnPath[r.id].append(p.id)
        
Rxn = defaultdict(list)

for x in model.reactions:
    Rxn[re.sub('[a-z]+','', x.reaction)].extend([[x.id,re.findall('[a-z]+', x.reaction)[0]]])
    
RxnID = int(re.findall('[0-9]+', sorted([x.id for x in model.reactions])[-1])[0])+1

ListOfMetFrom = {x.id: x.formula for x in model.metabolites}

Atom_ID = {'H2O': 'MAM02040', 'H': 'MAM02039', 'Fe': 'MAM01821', 'X': 'MAM01823', 'R': 'MAM03573', 'Na': 'MAM02519', 'K': 'MAM02200', 'Ca': 'MAM01413'}


rr = [x.id for x in model.reactions]           
m2 = [x.id for x in model.metabolites]
geneList = [x.id for x in model.genes]

print(f'{len(model.reactions)} reactions initially')
print(f'{len(model.metabolites)} metabolites initially')
print(f'{len(model.genes)} genes initially')

list_of_compartments = []
model2 = model.copy()


with model:
    for reac_count, x in enumerate(model.reactions[0:]): 

        x2 = model2.reactions.get_by_id(x.id)        

        print(f'Reaction({reac_count}/{len(rr)})')
        try:
            listOfgeneList5 = []
            if 'ec-code' in x2.annotation:
                EC = x.annotation['ec-code'] 
            else:
                EC = ''
            if type(EC) == str: 
                EC=[EC] 
            species = [s.id for s in x2.reactants] + [s.id for s in x2.products]
            species3 = [{s.split(' ')[0]: -1} if len(s.split(' ')) == 1 else {s.split(' ')[1]: -abs(float(s.split(' ')[0]))} for s in re.split(' --> | <=> ', re.sub('[a-z]+', '',x2.reaction))[0].split(' + ')] + [{s.split(' ')[0]: 1} if len(s.split(' ')) == 1 else {s.split(' ')[1]: float(s.split(' ')[0])} for s in re.split(' --> | <=> ', re.sub('[a-z]+', '',x2.reaction))[1].split(' + ')] 
            species3 = dict(ChainMap(*species3))
            mb = x2.check_mass_balance()
            bounds = x2.bounds
            gpr = x2.gpr
            annotation2 = x2.annotation
            locations = list(set([model2.metabolites.get_by_id(s).compartment for s in species]))
                      
            eq = x2.reaction
            eq_test = re.sub('<=>', '->', re.sub('-->','->',eq)).strip()
            eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()

            all_met_have_formula_test = min([1 if ListOfMetFrom[x] else 0 for x in species])
            if all_met_have_formula_test == 1:
                FromList=list()
                for y in species:
                    eq = re.sub(y, ListOfMetFrom[y], eq)
                    eq_test = re.sub(y, ListOfMetFrom[y], eq_test)

            
            if len(species) > 1 and all_met_have_formula_test == 1 and not test_reaction_balance(eq_test)[0]:  # first condition avoids exchange and sink reactions, second condition avoid to mass balance eq with some missing formula # Previous, worng: if len(FromList) > 1 and not '' in FromList: # MAR13082
                
                MB = mass_balance(eq, 'R')
                   
                NewSpecies = MB[6]+MB[7] # HH
                Stch = [-abs(x) for x in MB[0]]+MB[1]
                if MB[4] and MB[-1] != 2: 
                    [x2.add_metabolites({model.metabolites.get_by_id(Atom_ID[MB[4][0]]+locations[0]): MB[5][i]}) for i in range(len(MB[4]))]
                    species = [s.id for s in x2.reactants] + [s.id for s in x2.products] #
                if MB[-1] != 2: # Original: MB[-2]
                    for y in x2.metabolites:
                        x2.metabolites[y] = Stch[NewSpecies.index(y.formula)]
                        z = sympy.symbols('z')
                        expr = sympy.Eq(x2.metabolites[y]+z,Stch[NewSpecies.index(y.formula)])
                        sol = sympy.solve(expr)
                        x2.add_metabolites({model.metabolites.get_by_id(y.id): round(float(sol[0]),1)})

                    species3 = [{s.split(' ')[0]: -1} if len(s.split(' ')) == 1 else {s.split(' ')[1]: -abs(float(s.split(' ')[0]))} for s in re.split(' --> | <=> ', re.sub('[a-z]+', '',x2.reaction))[0].split(' + ')] + [{s.split(' ')[0]: 1} if len(s.split(' ')) == 1 else {s.split(' ')[1]: float(s.split(' ')[0])} for s in re.split(' --> | <=> ', re.sub('[a-z]+', '',x2.reaction))[1].split(' + ')] 
                    species3 = dict(ChainMap(*species3))
            
            if EC[0] and len(locations) == 1: # first: only reactions with EC can be analyzed, second: transport reactions are not analyzed
                
                for ec in EC:
                    print(ec)
                    if not ec in variables:
                        new_gpr = getGPR(ec, session) # new function 
                        if new_gpr[-1]:
                            new_locations = getLocation(new_gpr[3], new_gpr[1], new_gpr[2], 1, location_pkl_file, session) # original:(new_gpr[3], new_gpr[4], new_gpr[2], 1, location_pkl_file, session) # new function 
                            variables[ec].extend(new_gpr)
                            variables[ec].extend(new_locations)
                    else:
                        new_gpr = variables[ec][0:5]
                        print(new_gpr)
                        
                        new_locations = variables[ec][5:]
                        print(new_locations)
                    if new_gpr[-1]:
                        listOfgeneList5.append(variables[ec][5:])
                if listOfgeneList5 and not x2.reaction in RxnList:
                    print(listOfgeneList5)
                    geneList5 = meltGeneList(listOfgeneList5)
                    variables[x2.id]= geneList5
                    for CSL in [*variables[x2.id][0]]:                        
                        
                        CSL2 = CSL
                        
                        if not CSL in c: 
                            print()
                            print(0, CSL)
                            print()
                            CSL = 'cytosol'
                        
                        listAB = list(set(listB+listA))
                        if CSL in c:
                            ID = c[CSL]
                        elif len(re.sub(' $','',re.sub('^ ','',CSL)).split(" "))>1:
                            ID = (CSL.strip().split(" ")[0][0]+CSL.strip().split(" ")[1][0]).lower().replace(" ","")
                        else:
                            if len(CSL.split(" "))>1: 
                                ID = CSL[0:3].lower().replace(" ","")
                            else:
                                ID = CSL[0:2].lower().replace(" ","")
                        if not c.get(CSL) and ID in listAB:
                            r = re.compile(ID)
                            ID = ID+str(len(list(filter(r.match, listAB)))+1)
                        ComptoID2[CSL]=""
                        ComptoID2[CSL]+=ID
                        listA.append(ID)
                        listB.append(ID)
                        
                        if not ID in [rxn[1] for rxn in Rxn[re.sub('[a-z]+[0-9]*','', x2.reaction)]]:
                            reaction2 = Reaction('MAR'+str(RxnID+n))
                            reaction2.lower_bound = bounds[0]
                            reaction2.upper_bound = bounds[1]
                            for species2 in [[re.sub('[a-z][0-9]*',ComptoID2[CSL],x),x] for x in species]:
                                try:
                                    m = model2.metabolites.get_by_id(species2[0])
                                except KeyError:
                                    m = Metabolite(species2[0], charge=model2.metabolites.get_by_id(species2[1]).charge, formula=model2.metabolites.get_by_id(species2[1]).formula, name=model2.metabolites.get_by_id(species2[1]).name, compartment=ID)
                                    m.annotation = model2.metabolites.get_by_id(species2[1]).annotation

                                # Identify if the metabolite exists and if not add it to the model
                                if not m.id in model2.metabolites:
                                    model2.metabolites.add(m) # add metabolite m to the model

                                # Add compartment name if not exists
                                if not model2.compartments[ID]: 
                                    model2.compartments[ID] = CSL2

                                # Add the metabolite "species2" to the reaction x
                                reaction2.add_metabolites({m: species3[re.sub('[a-z]+[0-9]*', '',species2[0])]})
                 
                            if 'or' in variables[x2.id][2][CSL2] and 'and' in variables[x2.id][2][CSL2]:
                                reaction2.gene_reaction_rule = ' and '.join(['('+v+')' for v in variables[x2.id][2][CSL2].split(' and ')])
                            else:
                                reaction2.gene_reaction_rule = variables[x2.id][2][CSL2]
                                
                            reaction2.annotation = annotation2

                            model2.add_reactions([reaction2])

                            print(model2.reactions.get_by_id(reaction2.id))
                            print(len(model2.reactions))
                            for gene in re.split(' or | and ', reaction2.gene_reaction_rule):
                                
                                gene = gene.replace(')','').replace('(','')
                                if not gene in variables:
                                    ensemble = sorted(list(set(re.findall('ENS[A-Z][0-9]+', str(urllib.request.urlopen('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+gene).read())))))
                                    variables[gene] = ensemble
                                else: 
                                    ensemble = variables[gene]
                                    
                                model2.genes.get_by_id(gene).annotation['ensembl'] = ensemble
                                model2.genes.get_by_id(gene).annotation['hgnc.symbol'] = variables[x2.id][3][gene]

                            model2.groups.get_by_id(RxnPath[x2.id][0]).members.add(reaction2)

                            RxnList.append(reaction2.reaction)
                            n = n+1
                        else:
                            x2.gene_reaction_rule = variables[x2.id][2][CSL2]
                            for gene in re.split(' or | and ', x2.gene_reaction_rule):
                                gene = gene.replace(')','').replace('(','')
                                if not gene in variables:
                                    ensemble = sorted(list(set(re.findall('ENS[A-Z][0-9]+', str(urllib.request.urlopen('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+gene).read())))))
                                    variables[gene]=ensemble
                                else: 
                                    ensemble = variables[gene]

                                model2.genes.get_by_id(gene).annotation['hgnc.symbol']= variables[x2.id][3][gene]
                                model2.genes.get_by_id(gene).annotation['ensembl'] = ensemble
                        
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            print(x2.id)
            eList.append(x2.id)
            handle.write(x2.id+'\n')

        print('final',len(model.reactions))
        print('final2',len(model2.reactions))
      

    print()
    print(f'{len(model2.reactions)-len(rr)} new reactions')
    print([rxn.id for rxn in model2.reactions if not rxn.id in rr])
    
    print(f'{len(model2.metabolites)-len(m2)} new metabolites')
    print([m.id for m in model2.metabolites if not m.id in m2])
    
    print(f'{len(model2.genes)-len(geneList)} new genes')
    print([g.id for g in model2.genes if not g.id in geneList])
    
    print()
    print('An error was caused for the following reactions')
    
    print(eList)

    print('Saving variable environment')
    with open(variablesFile,'wb') as savedVariables:
        for variable in [variables]:
            dill.dump(variable,savedVariables,protocol=pickle.HIGHEST_PROTOCOL)
    print('Variable environment saved in '+variablesFile)

    print('Saving model')
    write_sbml_model(model, Output)
    write_sbml_model(model2, Output3)
   
    print('The model has been successfully written to '+Output+'\n')
    
    handle.close()


model3 = read_sbml_model(Output3)
model3_sec_copy = copy.deepcopy(model3)

########################### Evaluate and Correct the model ############################
model = read_sbml_model('model4_bb_2023-01-31.xml')

model4 = copy.deepcopy(model3)

## Check for duplicate reactions and eliminate

reactions_and_metabolites = [(x.id, str([(y.id, x.metabolites[y]) for y in x.metabolites])) for x in [x for x in model4.reactions]]

pattern = [str(x[1]) for x in reactions_and_metabolites]
pattern = [x for x in pattern if pattern.count(x) > 1]
pattern = list(set(pattern))

threshold = len(model.reactions)
duplicate_reactions = []
removed_reactions = []
for count, x in enumerate(pattern):
    print(count,'/',len(pattern))
    index = [i for i, y in enumerate([n[1] for n in reactions_and_metabolites]) if y == x]

    if len(index) > 1: # duplicate reactions

        duplicate_reactions.append(index)
        reactions_in_original_model = [x for x in index if x < threshold]
        
        if len(reactions_in_original_model) < 1: # if no reaction in the original model, then keep the reaction with the lowest index
            index = index[1:len(index)] 
        else: # otherwise, keep the reactions in the original model and delete the newly added reactions
            index = [x for x in index if x > threshold]

        for y in index:
            for g in model4.groups: # remove removed reactions from groups
                if reactions_and_metabolites[y][0] in [m.id for m in g.members]:
                    g.members.remove(model4.reactions.get_by_id(reactions_and_metabolites[y][0]))
            model4.reactions.remove(model4.reactions.get_by_id(reactions_and_metabolites[y][0]))
            removed_reactions.append(reactions_and_metabolites[y][0])


## Check for unconnected metabolites

metabolites_in_model = [x.id for x in model4.metabolites]
metabolites_in_reactions = [[y.id for y in x.metabolites.keys()] for x in model4.reactions]

metabolites_in_reactions = []
for x in model4.reactions:
    for y in x.metabolites:
        metabolites_in_reactions.append(y.id)
metabolites_in_reactions = list(set(metabolites_in_reactions))

isolated_metabolites = [x for x in metabolites_in_model if not x in metabolites_in_reactions]
for x in isolated_metabolites:
    model4.metabolites.remove(x)


## Check for umbalanced reactions

list_reaction_balance = []
for x in model4.reactions:
    print(x.id)
    eq = x.reaction
    eq_test = re.sub('<=>', '->', re.sub('-->','->',eq)).strip()
    eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()
    species = [s.id for s in x.reactants] + [s.id for s in x.products]    
    try:
        all_met_have_formula_test = min([1 if ListOfMetFrom[x] else 0 for x in species])
    except:
        all_met_have_formula_test = 0
    if all_met_have_formula_test == 1:
        FromList = list()
        for y in species:
           eq_test = re.sub(y, ListOfMetFrom[y], eq_test)
           eq = re.sub(y, ListOfMetFrom[y], eq)           
        if 1 < len(species) < 25:   
            mass_balance_test = test_reaction_balance(eq_test)   
            if not mass_balance_test[0]:
                MB = mass_balance(eq, 'R')
                if model.reactions.get_by_id(x.id):
                    original_model_reaction = model.reactions.get_by_id(x.id).reaction
                else:
                    original_model_reaction = ''
                list_reaction_balance.append([x.id, mass_balance_test, MB, original_model_reaction])


# Writte the sbml model:

write_sbml_model(model4, Output4) #Human1.5


########################### Generate THG_beta1 from THG_beta2 ############################

model5 = copy.deepcopy(model4)

# Reactions

list_of_reactions_in_thgb2_not_in_thgb1 = [x.id for x in model4.reactions if not x.id in [y.id for y in model.reactions]]

for x in list_of_reactions_in_thgb2_not_in_thgb1:
    for g in model5.groups: # remove removed reactions from groups
        if x in [m.id for m in g.members]:
            g.members.remove(model5.reactions.get_by_id(x))
    model5.reactions.remove(x)

# Metabolites

list_of_metabolites_in_thgb2_not_in_thgb1 = [x.id for x in model4.metabolites if not x.id in [y.id for y in model.metabolites]]

for x in list_of_metabolites_in_thgb2_not_in_thgb1: model5.metabolites.remove(x)

## Check for umbalanced reactions

list_reaction_balance = []
for x in model5.reactions:
    print(x.id)
    eq = x.reaction
    eq_test = re.sub('<=>', '->', re.sub('-->','->',eq)).strip()
    eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()
    species = [s.id for s in x.reactants] + [s.id for s in x.products]    
    try:
        all_met_have_formula_test = min([1 if ListOfMetFrom[x] else 0 for x in species])
    except:
        all_met_have_formula_test = 0
    if all_met_have_formula_test == 1:
        FromList = list()
        for y in species:
           eq_test = re.sub(y, ListOfMetFrom[y], eq_test)
           eq = re.sub(y, ListOfMetFrom[y], eq)           
        if 1 < len(species) < 25:   
            mass_balance_test = test_reaction_balance(eq_test)   
            if not mass_balance_test[0]:
                MB = mass_balance(eq, 'R')
                if model.reactions.get_by_id(x.id):
                    original_model_reaction = model.reactions.get_by_id(x.id).reaction
                else:
                    original_model_reaction = ''
                list_reaction_balance.append([x.id, mass_balance_test, MB, original_model_reaction])

# Writte the sbml model:

write_sbml_model(model5, Output5) #Human1.2

# Test if sbml models can be read

test1 = read_sbml_model(Output4)
test2 = read_sbml_model(Output5)
