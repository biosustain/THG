from cobra.io import read_sbml_model,write_sbml_model
import os
import sys
from pathlib import Path
# -*- coding: utf-8 -*-
import urllib.request, urllib.error, urllib.parse#,cookielib
import re
import urllib.request, urllib.parse, urllib.error
import requests
#import numpy as np
#from scipy import linalg
import copy
import time
#from itertools import groupby
import traceback
import itertools
import pubchempy as pcp
import string
import pickle 
from collections import defaultdict
from itertools import zip_longest
#from Class import * ##### for 
#from Pattern import * ###### for 
#from Function import * ###### for 
from equations_build_model import *
from cobra import Model, Reaction, Metabolite
#from Programa_3 import *
from collections import ChainMap
import dill
from function_build_model import *


from gpr.auth_gpr import getGPR, setup_biocyc_session
# 👇 import the addtional function
from gpr.getLocation import getLocation
session = setup_biocyc_session()
if len(sys.argv) > 1:
    location_pkl_file = sys.argv[1]
else:
    location_pkl_file = 'pkl/Human-variables.pkl'
    print(f"No location file provided, using '{location_pkl_file}'")



model = read_sbml_model('models/H1.1.xml')
model2 = model


variablesFile = 'files/variablesFile.pkl' 
Output = f'models/Output_{Path(location_pkl_file).stem}.xml'  
print(f"Will write model to {Output}")
Output2 = 'files/Rxn2Fix'

#if os.path.exists(variablesFile):
#    with open(variablesFile, 'rb') as handle:
#        variables = pickle.load(handle)
#else:
#    variables = defaultdict()
    
handle = open(Output2, 'w')



#with open("pkl/Compvariable", 'rb') as handle: ### the most recent dict from the database 
#    c = pickle.load(handle)

c= compartment_file_to_dict() 

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

#idx = [x.id for x in model.reactions].index()
           
rr = [x.id for x in model.reactions]           
m2 = [x.id for x in model.metabolites]
geneList = [x.id for x in model.genes]

print(f'{len(model.reactions)} reactions initially')
print(f'{len(model.metabolites)} metabolites initially')
print(f'{len(model.genes)} genes initially')

#idx = [x.id for x in model.reactions].index('MAR07161') # IM to HH: What is this reaction???
with model:
    for reac_count, x in enumerate(model2.reactions):   
        print(f'Reaction({reac_count+1}/{len(rr)})')
        try:
            listOfgeneList5 = []
            if 'ec-code' in x.annotation:
                EC = x.annotation['ec-code'] 
            else:
                EC = ''
            if type(EC) == str: EC=[EC] 
            species = [x.id for x in x.reactants] + [x.id for x in x.products]
            species3 = [{x.split(' ')[0]: -1} if len(x.split(' ')) == 1 else {x.split(' ')[1]: -abs(float(x.split(' ')[0]))} for x in re.split(' --> | <=> ', re.sub('[a-z]+', '',x.reaction))[0].split(' + ')] + [{x.split(' ')[0]: 1} if len(x.split(' ')) == 1 else {x.split(' ')[1]: float(x.split(' ')[0])} for x in re.split(' --> | <=> ', re.sub('[a-z]+', '',x.reaction))[1].split(' + ')] 
            species3 = dict(ChainMap(*species3))
            mb = x.check_mass_balance()
            bounds = x.bounds
            gpr = x.gpr
            annotation2 = x.annotation
            locations = list(set([model.metabolites.get_by_id(x).compartment for x in species]))
                      
            eq = x.reaction  
            eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()
            FromList=list()
            for y in species:
                eq = re.sub(y, ListOfMetFrom[y], eq)

            all_met_have_formula_test = min([1 if x in ListOfMetFrom.keys() else 0 for x in species]) # 0:no, 1:yes
            
            if len(species) > 1 and all_met_have_formula_test == 1: # first condition avoids exchange and sink reactions, second condition avoid to mass balance eq with some missing formula # Previous, worng: if len(FromList) > 1 and not '' in FromList: # MAR13082
                
                MB = RxnBalance2(eq, 'R')
                   
                #NewSpecies = MB[4]+MB[5] # Original: MB[6]+MB[7] # IM
                NewSpecies = MB[6]+MB[7] # HH
                Stch = [-abs(x) for x in MB[0]]+MB[1]
                if MB[4] and MB[-1] != 2: # Original: MB[-2]RxnBalance2 output: MB[0], MB[1], AddH, AddH2O, NewSpecies[0], NewSpecies[1], AllSpecies[0], AllSpecies[1], eq, eq_init, TestOfBalance
                    [x.add_metabolites({model.metabolites.get_by_id(Atom_ID[MB[4][0]]+locations[0]): MB[5][i]}) for i in range(len(MB[4]))]
                    species = [x.id for x in x.reactants] + [x.id for x in x.products] #
                if MB[-1] != 2: # Original: MB[-2]
                    for y in x.metabolites:
                        x.metabolites[y] = Stch[NewSpecies.index(y.formula)]
                        z = sympy.symbols('z')
                        expr = sympy.Eq(x.metabolites[y]+z,Stch[NewSpecies.index(y.formula)])
                        sol = sympy.solve(expr)
                        x.add_metabolites({model.metabolites.get_by_id(y.id): round(float(sol[0]),1)})
                    species3 = [{x.split(' ')[0]: -1} if len(x.split(' ')) == 1 else {x.split(' ')[1]: -abs(float(x.split(' ')[0]))} for x in re.split(' --> | <=> ', re.sub('[a-z]+', '',x.reaction))[0].split(' + ')] + [{x.split(' ')[0]: 1} if len(x.split(' ')) == 1 else {x.split(' ')[1]: float(x.split(' ')[0])} for x in re.split(' --> | <=> ', re.sub('[a-z]+', '',x.reaction))[1].split(' + ')] 
                    species3 = dict(ChainMap(*species3))
            
            
            if EC[0] and len(locations) == 1: # first: only reactions with EC can be analyzed, second: transport reactions are not analyzed
                
                for ec in EC:
                    print(ec)
                    if not ec in variables:
                        #new_gpr = getGPR22(ec,20) # old function
                        new_gpr = getGPR(ec, session) # new function 
                        if new_gpr[-1]:
                            #print(new_gpr)
                            #new_locations = getLocation_old(new_gpr[3], new_gpr[4], new_gpr[2], 20) # old function 
                            new_locations = getLocation(new_gpr[3], new_gpr[1], new_gpr[2], 1, location_pkl_file, session) # original:(new_gpr[3], new_gpr[4], new_gpr[2], 1, location_pkl_file, session) # new function 
                            variables[ec].extend(new_gpr)
                            variables[ec].extend(new_locations)
                            #print(new_locations)
                    else:
                        new_gpr = variables[ec][0:5]
                        print(new_gpr)
                        
                        new_locations = variables[ec][5:]
                        print(new_locations)
                    if new_gpr[-1]:
                        listOfgeneList5.append(variables[ec][5:])
                if listOfgeneList5 and not x.reaction in RxnList:
                    print(listOfgeneList5)
                    geneList5 = meltGeneList(listOfgeneList5)
                    variables[x.id]= geneList5
                    for CSL in [*variables[x.id][0]]:
                        
                        
                        CSL2 = CSL
                        CSL = CSL[0].upper()+CSL[1:]
                        
                        #if not CSL in c: 
                        #    print()
                        #    print(0, CSL)
                        #    print()
                        #    CSL = 'Cytosol'
                        
                        listAB = list(set(listB+listA))
                        if CSL in c:
                            ID = c[CSL]
                        elif len(re.sub(' $','',re.sub('^ ','',CSL)).split(" "))>1: # would be equivalent to len(CSL.strip().split(" "))>1
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
                        
                        
                        
                        
                        if not ID in [x[1] for x in Rxn[re.sub('[a-z]+[0-9]*','', x.reaction)]]:
                            reaction2 = Reaction('MAR'+str(RxnID+n))
                            reaction2.lower_bound = bounds[0]
                            reaction2.upper_bound = bounds[1]
                            for species2 in [[re.sub('[a-z][0-9]*',ComptoID2[CSL],x),x] for x in species]:
                                try:
                                    m = model.metabolites.get_by_id(species2[0])
                                    reaction2.add_metabolites({m: species3[re.sub('[a-z]+[0-9]*', '',species2[0])]})
                                except KeyError:
                                    m=Metabolite(species2[0], charge=model.metabolites.get_by_id(species2[1]).charge, formula=model.metabolites.get_by_id(species2[1]).formula, name=model.metabolites.get_by_id(species2[1]).name, compartment=model.metabolites.get_by_id(species2[1]).compartment)
                                    m.annotation = model.metabolites.get_by_id(species2[1]).annotation
                                    reaction2.add_metabolites({m: species3[re.sub('[a-z]+[0-9]*', '',species2[0])]})
                                    
                            if 'or' in variables[x.id][2][CSL2] and 'and' in variables[x.id][2][CSL2]:
                                reaction2.gene_reaction_rule = ' and '.join(['('+x+')' for x in variables[x.id][2][CSL2].split(' and ')])
                            else:
                                reaction2.gene_reaction_rule = variables[x.id][2][CSL2]
                                
                            reaction2.annotation = annotation2
                            model.add_reactions([reaction2])
                            
                            for gene in re.split(' or | and ', reaction2.gene_reaction_rule):
                                
                                gene = gene.replace(')','').replace('(','')
                                if not gene in variables:
                                    ensemble = sorted(list(set(re.findall('ENS[A-Z][0-9]+', str(urllib.request.urlopen('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+gene).read())))))
                                    variables[gene]=ensemble
                                else: 
                                    ensemble = variables[gene]
                                    
                                model.genes.get_by_id(gene).annotation['ensembl'] = ensemble
                                model.genes.get_by_id(gene).annotation['hgnc.symbol']=variables[x.id][3][gene]
                                
                            model.groups.get_by_id(RxnPath[x.id][0]).members.add(reaction2)
                            RxnList.append(reaction2.reaction)
                            n = n+1
                        else:
                            x.gene_reaction_rule = variables[x.id][2][CSL2]
                            for gene in re.split(' or | and ', x.gene_reaction_rule):
                                gene = gene.replace(')','').replace('(','')
                                if not gene in variables:
                                    ensemble = sorted(list(set(re.findall('ENS[A-Z][0-9]+', str(urllib.request.urlopen('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+gene).read())))))
                                    variables[gene]=ensemble
                                else: 
                                    ensemble = variables[gene]
                                model.genes.get_by_id(gene).annotation['hgnc.symbol']= variables[x.id][3][gene]
                                model.genes.get_by_id(gene).annotation['ensembl'] = ensemble

                        
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            print(x.id)
            eList.append(x.id)
            handle.write(x.id+'\n')
    print()
    print(f'{len(model.reactions)-len(rr)} new reactions')
    print([x.id for x in model.reactions if not x.id in rr])
    
    print(f'{len(model.metabolites)-len(m2)} new metabolites')
    print([x.id for x in model.metabolites if not x.id in m2])
    
    print(f'{len(model.genes)-len(geneList)} new genes')
    print([x.id for x in model.genes if not x.id in geneList])
    
    print()
    print('An error was caused for the following reactions')
    
    print(eList)
    #print()
    print('Saving variable environment')
    with open(variablesFile,'wb') as savedVariables:
        for variable in [variables]:
            dill.dump(variable,savedVariables,protocol=pickle.HIGHEST_PROTOCOL)
    print('Variable environment saved in '+variablesFile)
    #print()

    print('Saving model')
    write_sbml_model(model, Output)
    print('The model has been successfully written to '+Output+'\n')
    
    handle.close()

