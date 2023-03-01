# -*- coding: utf-8 -*-

# Libraries
import re
import pickle
import gzip, pickle
from ast import And, BoolOp, Expression, Name, Or
from functools import reduce
from typing import List
import cobra
from cobra.io import load_model
from cobra.io import validate_sbml_model
import cobra.io
import cobra.core
import cobra.core.metabolite
from cobra import Model, Reaction, Metabolite
from cobra.core.gene import GPR, ast_parse
import pandas as pd
import copy
from math import isnan
from functools import reduce
import numpy as np
import difflib
from cobra.io import read_sbml_model,write_sbml_model
import memote
from collections import defaultdict
import memote.support.consistency as consistency
import memote.support.consistency_helpers as con_helpers
from memote.utils import annotate, get_ids, truncate, wrapper
import logging
import memote.support.consistency_helpers as con_helpers
import memote.support.helpers as helpers
import pytest
from functions_network_consistency import *
from functions_mass_balance import *
from collections import ChainMap


# Functions
# Used in several steps
'Calculate Jaccard Index'
def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union


def test_reaction_balance(eq):  
	eq=" " + eq # add an extra space at the beginning to avoid problems when defining StCoeff of the first metabolite using reduce func    
	Ls = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
	#Mass Balance
	Ss,Es,i,ii=defaultdict(list),[],1,defaultdict(list)
	for p in eq.split('->'):
		for kk in p.split('+'):
			StCoeff =[reduce(lambda i, j: i + j, kk[0 : re.search("[A-Z]",kk).start()])] # This line improves the previous AtomCount function allowing to account for the stoichiometric coeff
			StCoeff = 1 if StCoeff[0] == ' ' else float(StCoeff[0]) #Stoichiometric Coeff improvement
			c = [Ls.pop(0), 1]
			for e,m in re.findall('([A-Z][a-z]?)([0-9]*)',kk):              
				m = StCoeff if m == '' else int(m)*StCoeff #Stoichiometric Coeff improvement
				d = [c[0],c[1]*m*i]
				Ss[e][:0],Es[:0] = [d],[[e,d]]            
				ii[e].append(d[1]) 
		i=-1  
	MassBalance = sum([abs(round(sum(x),2)) for x in ii.values()])
	if MassBalance == 0:   
		Species = [x.strip() for x in [x for x in eq.split('->')][0].split('+')],[x.strip() for x in [x for x in eq.split('->')][1].split('+')]
		i=0
		for x in Species[0]:
			if re.findall('[A-Z]',x[0]): Species[0][i] = '1 '+ x
			else: Species[0][i] = x
			i += 1
		i=0
		for x in Species[1]:
			if re.findall('[A-Z]',x[0]): Species[1][i] = '1 '+ x
			else: Species[1][i] = x
			i += 1
		SubsStc = [x.split(' ')[0] for x in [x.strip() for x in Species[0]]]
		ProdStc = [x.split(' ')[0] for x in [x.strip() for x in Species[1]]]
		SubsStch = [float(x) for x in SubsStc]
		ProdStch = [float(x) for x in ProdStc]            
	else:
		SubsStch = ''
		ProdStch = ''   
	return SubsStch , ProdStch, eq


def imbalance_test (x, ListOfMetFrom):
    '''
    Input:
        x = reaction from model
        ListOfMetFrom = {x.id: x.formula for x in model.metabolites}
    
    '''
    species = [x.id for x in x.reactants] + [x.id for x in x.products]
    all_met_have_formula_test = min([1 if ListOfMetFrom[x] else 0 for x in species]) # 0:no, 1:yes
    if 1 < len(species) < 20 and all_met_have_formula_test == 1:
        eq = x.reaction
        eq = re.sub('<=>', '->', re.sub('-->','->',eq)).strip()
        for y in species:
            eq = re.sub(y, ListOfMetFrom[y], eq)
        mb = test_reaction_balance(eq)
        return mb


# Step 1: transfer metabolites from source network to target network
'Compare two lists of metabolites from two different reactions'
def compare_metabolite_groups (list_1,list_2,model_1,model_2):  # Takes into account name, formula and the common attributes defined in annotation
    """
    For two sets of metabolites, the function takes the common atributes and compare them using jJaccard index metrics to determine if the two sets are equal or not

    Inputs
    ----------
    list1 and list2 : list of metabolites IDs (ID used to be identified in the respective models) (i.e.: ['MAM01306m', 'MAM01370m'])
    model_1 and model_2: metabolic models in cobra format corresponding to list_1 and list_2 respectively
    
    Output
    -------
    overlap: metric to evaluate the overlap between the two lists based on Jaccard Index. 
        0: no overlap
        1: list_1 and list_2 are identical
    item_comparison: list scoring similarlity between all the items in both lists
        0: 0% similar
        1: 100% similar
    """
    if isinstance(list_1,str): list_1 = [list_1]
    if isinstance(list_2,str): list_2 = [list_2]
    metabolite_list_1 = [model_1.metabolites.get_by_id(t) for t in list_1]
    metabolite_list_2 = [model_2.metabolites.get_by_id(t) for t in list_2]
    # compare all the elements in list_1 against all the elements in list_2 using string comparison (name and formula)
    if [1 for x in metabolite_list_1+metabolite_list_2 if x.name == None]:
        check_name = [0]*(len(list_1)*len(list_2))
    else:
        check_name = [difflib.SequenceMatcher(None,y.lower(),x.lower()).ratio() for y in [model_1.metabolites.get_by_id(z).name for z in list_1] for x in [model_2.metabolites.get_by_id(t).name for t in list_2]]   
    if [1 for x in metabolite_list_1+metabolite_list_2 if x.formula == None]:
        check_formula = [0]*(len(list_1)*len(list_2))
    else:
        check_formula = [difflib.SequenceMatcher(None,y,x).ratio() for y in [model_1.metabolites.get_by_id(z).formula for z in list_1] for x in [model_2.metabolites.get_by_id(t).formula for t in list_2]]
    item_comparison = max(check_name,check_formula)
    overlap = min(item_comparison)
    # try to find a match by using the metabolites attributes (IDs)
    r = list(reduce(lambda i, j: i & j, (set(x.annotation.keys()) for x in metabolite_list_1+metabolite_list_2))) # list of common attributes that will be used to make the comparison
    for x in r: 
        if x != 'sbo': # sbo id is common for all the metabolites
            xth_list_1 = [z.annotation.get(x) for z in metabolite_list_1]
            xth_list_2 = [z.annotation.get(x) for z in metabolite_list_2]
            for p in xth_list_1:
                for q in xth_list_2:
                    overlap_test = jaccard(p, q)
                    if overlap_test > overlap: 
                        overlap = overlap_test
                        item_comparison = [difflib.SequenceMatcher(None,y,x).ratio() for y in xth_list_1 for x in xth_list_2]
    return overlap, item_comparison

'Enrich metabolite list and annotation in network_1 from network_2'
def network_metabolites_merge (network_1, network_2): # target network, source network
    """
    Metabolites from two metabolic models are compared based on name, formula and other IDs. If a match is found the information from network_2 replaces the infromation in network_1

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated metabolite information from network_2 for the overlaped metabolites and new metabolites from network_2 in case they were not annotated in network_1 previously
    n1_n2_equivalent_metabolites : dictionary of metabolite IDs equivalent between network_1 and network_2. key: metabolite ID in network_2, value: metabolite ID in network_1
    """
    metabolites_in_network_1 = [x.id for x in network_1.metabolites]
    metabolites_in_network_2 = [x.id for x in network_2.metabolites]
    unique_metabolites_in_network_1 = set([re.sub(x.compartment+'$', '', x.id) for x in network_1.metabolites])
    unique_metabolites_in_network_2 = set([re.sub(x.compartment+'$', '', x.id) for x in network_2.metabolites])
    unique_pairs_to_test = [[x,y] for x in unique_metabolites_in_network_1 for y in unique_metabolites_in_network_2]
    n1_n2_equivalent_metabolites = {} # dictionary created here and used when upgrading reactions, the metabolite ids in the source reactions (from n2) will be changed to fit with the annotation in the target network (n1)
    for x in network_1.metabolites: # x metabolite in target network
        unique_x = re.sub(x.compartment+'$', '', x.id) # metabolite id with no compartment ident
        for y in network_2.metabolites: # y metabolite in source network
            unique_y = re.sub(y.compartment+'$', '', y.id) # metabolite id with no compartment ident
            pair_2_string = [unique_x, unique_y] # pair of unique metabolites to compare
            if pair_2_string in unique_pairs_to_test: # check if the pair is in the list to compare, if not skip (already compared)
                print(1,x,y)
                unique_x_y_pair = unique_pairs_to_test.pop(unique_pairs_to_test.index(pair_2_string)) # remove the pair to compare from the list of pairs to compare to not repeat the process      
                x_in_n1 = [x for x in metabolites_in_network_1 if re.search(unique_x, x)] # unique metabolite in all the compartments in n1
                y_in_n2 = [x for x in metabolites_in_network_2 if re.search(unique_y, x)] # unique metabolite in all the compartments in n2
                ref_target_metabolite =  network_1.metabolites.get_by_id(x_in_n1[0])
                ref_source_metabolite =  network_2.metabolites.get_by_id(y_in_n2[0])
                source_vs_target = compare_metabolite_groups(x.id,y.id,network_1,network_2)[0]
                if source_vs_target == 1 and unique_x in unique_metabolites_in_network_1 and unique_y in unique_metabolites_in_network_2:
                    print(2,x,y)
                    n1_n2_equivalent_metabolites[unique_y] = unique_x # {metabolite_n2 : metabolite_n1} list of equivalent metabolite id between network_1 and network_2
                    for x in x_in_n1: # replace the information of metabolites in target network (n1) by the right information in source network (n2)
                        network_1.metabolites.get_by_id(x).annotation = copy.deepcopy(ref_source_metabolite.annotation)
                        network_1.metabolites.get_by_id(x).charge = copy.deepcopy(ref_source_metabolite.charge)
                        network_1.metabolites.get_by_id(x).elements = copy.deepcopy(ref_source_metabolite.elements)
                        network_1.metabolites.get_by_id(x).formula = copy.deepcopy(ref_source_metabolite.formula)
                        network_1.metabolites.get_by_id(x).name = copy.deepcopy(ref_source_metabolite.name)
                        network_1.metabolites.get_by_id(x).notes = copy.deepcopy(ref_source_metabolite.notes)
                    metabolite_compartment_n1 = [x.compartment for x in [network_1.metabolites.get_by_id(x) for x in x_in_n1]]
                    for y in y_in_n2: # add metabolites in compartments not annotated in the target network 
                        yth_metabolite = network_2.metabolites.get_by_id(y)
                        if not yth_metabolite.compartment in metabolite_compartment_n1: 
                            print(3,x,y)
                            network_1.add_metabolites(yth_metabolite)# if the yth compartment/specific metabolite in the source network is not in the target network, add it 
                            network_1.metabolites.get_by_id(y).id = unique_x+yth_metabolite.compartment # change the id to be consistent with the rest                        
                    unique_metabolites_in_network_1.remove(unique_x)
                    unique_metabolites_in_network_2.remove(unique_y)
    for x in unique_metabolites_in_network_2: # add to target network (n1) all the metabolites in the source network (n2) with no equivalent metabolite in n1
        x = list(unique_metabolites_in_network_2)[0]
        metabolite_to_add = [y for y in metabolites_in_network_2 if re.search(x, y)]
        for y in metabolite_to_add: # add metabolites from n2 to n1
            yth_metabolite = network_2.metabolites.get_by_id(y)
            network_1.add_metabolites(yth_metabolite)
    return network_1, n1_n2_equivalent_metabolites


def compare_metabolite_groups_2 (list_1,list_2,model_1,model_2): # Only takes into account name and formula
    """
    For two sets of metabolites, the function takes the common atributes and compare them using jJaccard index metrics to determine if the two sets are equal or not

    Inputs
    ----------
    list1 and list2 : list of metabolites IDs (ID used to be identified in the respective models) (i.e.: ['MAM01306m', 'MAM01370m'])
    model_1 and model_2: metabolic models in cobra format corresponding to list_1 and list_2 respectively
    
    Output
    -------
    overlap: metric to evaluate the overlap between the two lists based on Jaccard Index. 
        0: no overlap
        1: list_1 and list_2 are identical
    item_comparison: list scoring similarlity between all the items in both lists
        0: 0% similar
        1: 100% similar
    """
    if isinstance(list_1,str): list_1 = [list_1]
    if isinstance(list_2,str): list_2 = [list_2]
    metabolite_list_1 = [model_1.metabolites.get_by_id(t) for t in list_1]
    metabolite_list_2 = [model_2.metabolites.get_by_id(t) for t in list_2]
    # compare all the elements in list_1 against all the elements in list_2 using string comparison (name and formula)
    if [1 for x in metabolite_list_1+metabolite_list_2 if x.name == None]:
        check_name = [0]*(len(list_1)*len(list_2))
    else:
        check_name = [difflib.SequenceMatcher(None,y.lower(),x.lower()).ratio() for y in [model_1.metabolites.get_by_id(z).name for z in list_1] for x in [model_2.metabolites.get_by_id(t).name for t in list_2]]   
    if [1 for x in metabolite_list_1+metabolite_list_2 if x.formula == None]:
        check_formula = [0]*(len(list_1)*len(list_2))
    else:
        check_formula = [difflib.SequenceMatcher(None,y,x).ratio() for y in [model_1.metabolites.get_by_id(z).formula for z in list_1] for x in [model_2.metabolites.get_by_id(t).formula for t in list_2]]
    overlap = max(check_name,check_formula)
    return overlap


def network_metabolites_merge_2 (network_1, network_2): # target network, source network # This is the most efficient, however it is deprecated for the moment as it gives some false matches due to some miss annotation in network 2. Once it is corrected in the DB building, this function will be reactivated
    """
    Metabolites from two metabolic models are compared based on name, formula and other IDs. If a match is found the information from network_2 replaces the infromation in network_1

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated metabolite information from network_2 for the overlaped metabolites and new metabolites from network_2 in case they were not annotated in network_1 previously
    n1_n2_equivalent_metabolites : dictionary of metabolite IDs equivalent between network_1 and network_2. key: metabolite ID in network_2, value: metabolite ID in network_1
    """
    # Attributes
    def attributes (network):
        list_of_metabolite_attributes = [] # List of all the metabolite attributes in network_1
        for s in network.metabolites:
            list_of_metabolite_attributes = list_of_metabolite_attributes + [x for x in s.annotation]
        list_of_metabolite_attributes = list(set(list_of_metabolite_attributes)) 
        attributes = {} # dictionary containing all the attributes of all the metabolites in network_1
        for attribute in list_of_metabolite_attributes:
            attributes[attribute] = {}
            for s in network.metabolites:
                if attribute in s.annotation:
                    attributes[attribute][str(s.annotation.get(attribute))] = re.sub(str(s.compartment)+"$","",s.id)
        return attributes

    list_of_metabolite_attributes_n1 = attributes (network_1) # List of all the metabolite attributes in network_1
    list_of_metabolite_attributes_n2 = attributes (network_2) # List of all the metabolite attributes in network_2
    list_of_common_attributes = list(list_of_metabolite_attributes_n1.keys() & list_of_metabolite_attributes_n2.keys())

    # Unique metabolites and compartments
    def unique_metabolites_non_unique_metabolites_and_compartments (network):
        unique_metabolites_and_compartments = {}
        for m in network.metabolites:
            mth_unique_metabolite = re.sub(str(m.compartment)+"$","",m.id)
            if not mth_unique_metabolite in unique_metabolites_and_compartments.keys():
                unique_metabolites_and_compartments[mth_unique_metabolite] = {} #{'compartment', 'non_unique_metabolites'}
                unique_metabolites_and_compartments[mth_unique_metabolite]['compartment'] = []
                unique_metabolites_and_compartments[mth_unique_metabolite]['non_unique_metabolites'] = []
            unique_metabolites_and_compartments[mth_unique_metabolite]['compartment'].append(m.compartment)
            unique_metabolites_and_compartments[mth_unique_metabolite]['non_unique_metabolites'].append(m.id)
        return unique_metabolites_and_compartments

    unique_metabolites_and_compartments_n1 = unique_metabolites_non_unique_metabolites_and_compartments(network_1)
    unique_metabolites_and_compartments_n2 = unique_metabolites_non_unique_metabolites_and_compartments(network_2)

    # Equivalent metabolites between n1 and n2
    n1_n2_equivalent_unique_metabolites = [] # dictionary created here and used when upgrading reactions, the metabolite ids in the source reactions (from n2) will be changed to fit with the annotation in the target network (n1)
    for x in list_of_common_attributes: # first try with the attributes
        if x != 'sbo': # avoid sbo attribute because it is the same for all the metabolites
            xth_n1_n2_intersection = list(list_of_metabolite_attributes_n1[x].keys() & list_of_metabolite_attributes_n2[x].keys())
            for y in xth_n1_n2_intersection:
                n1_n2_equivalent_unique_metabolites.append([list_of_metabolite_attributes_n1[x][y],list_of_metabolite_attributes_n2[x][y]])

    remaining_metabolites_n1 = [x for x in list(unique_metabolites_and_compartments_n1.keys()) if x not in list(set([x[0] for x in n1_n2_equivalent_unique_metabolites]))]
    remaining_metabolites_n2 = [x for x in list(unique_metabolites_and_compartments_n2.keys()) if x not in list(set([x[1] for x in n1_n2_equivalent_unique_metabolites]))]

    for m2 in remaining_metabolites_n2:
        m2th_non_unique_id = unique_metabolites_and_compartments_n2[m2]['non_unique_metabolites'][0]
        for m1 in remaining_metabolites_n1:
            m1th_non_unique_id = unique_metabolites_and_compartments_n1[m1]['non_unique_metabolites'][0]
            # Compare m2th and m1th
            m2th_vs_m1th = compare_metabolite_groups_2 (m1th_non_unique_id,m2th_non_unique_id,network_1,network_2)
            if m2th_vs_m1th == 1:
                print(m1,m2)
                n1_n2_equivalent_unique_metabolites.append(m1,m2)
                remaining_metabolites_n1 = [x for x in remaining_metabolites_n1 if x != m1]
                break

    # Enrich network 1 with the equivalent metabolites in network 2
    for x in n1_n2_equivalent_unique_metabolites: #enrich overlaped metabolites in network_1
        try:
            source_metabolite = network_2.metabolites.get_by_id(unique_metabolites_and_compartments_n2[x[1]]['non_unique_metabolites'][0])
            new_compartments = [c for c in unique_metabolites_and_compartments_n2[x[1]]['compartment'] if not c in unique_metabolites_and_compartments_n1[x[0]]['compartment']]      
            
            for y in unique_metabolites_and_compartments_n1[x[0]]['non_unique_metabolites']: # update the overlaped metabolites
                for a in network_1.metabolites.get_by_id(y).annotation.keys(): # annotation
                    if a in source_metabolite.annotation.keys():
                        network_1.metabolites.get_by_id(y).annotation[a] = source_metabolite.annotation[a]
                network_1.metabolites.get_by_id(y).charge = copy.deepcopy(source_metabolite.charge)
                network_1.metabolites.get_by_id(y).elements = copy.deepcopy(source_metabolite.elements)
                network_1.metabolites.get_by_id(y).formula = copy.deepcopy(source_metabolite.formula)
                network_1.metabolites.get_by_id(y).name = copy.deepcopy(source_metabolite.name)
                network_1.metabolites.get_by_id(y).notes = copy.deepcopy(source_metabolite.notes)

            template_metabolite = copy.deepcopy(network_1.metabolites.get_by_id(y))

            for z in new_compartments: # add metabolites in compartments not annotated in the target network 
                zth_metabolite = copy.deepcopy(template_metabolite)
                zth_metabolite.id = re.sub(zth_metabolite.compartment+'$', z, zth_metabolite.id)
                zth_metabolite.compartment = z 
                network_1.add_metabolites(zth_metabolite)
        except:
            continue

    # Add the non-overlapped metabolites
    remaining_metabolites_n2 = [x for x in list(unique_metabolites_and_compartments_n2.keys()) if x not in list(set([x[1] for x in n1_n2_equivalent_unique_metabolites]))]

    for x in remaining_metabolites_n2:
        non_unique_metabolites = unique_metabolites_and_compartments_n2[x]['non_unique_metabolites']
        for y in non_unique_metabolites:
            xyth_metabolite = copy.deepcopy(network_2.metabolites.get_by_id(y))
            network_1.add_metabolites(xyth_metabolite)

    return network_1, n1_n2_equivalent_unique_metabolites


def network_metabolites_merge_3 (network_1, network_2): # target network, source network
    """
    Metabolites from two metabolic models are compared based on name, formula and other IDs. If a match is found the information from network_2 replaces the infromation in network_1

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated metabolite information from network_2 for the overlaped metabolites and new metabolites from network_2 in case they were not annotated in network_1 previously
    n1_n2_equivalent_metabolites : dictionary of metabolite IDs equivalent between network_1 and network_2. key: metabolite ID in network_2, value: metabolite ID in network_1
    """
    # Attributes
    def attributes (network):
        list_of_metabolite_attributes = [] # List of all the metabolite attributes in network_1
        for s in network.metabolites:
            list_of_metabolite_attributes = list_of_metabolite_attributes + [x for x in s.annotation]
        list_of_metabolite_attributes = list(set(list_of_metabolite_attributes)) 
        attributes = {} # dictionary containing all the attributes of all the metabolites in network_1
        for attribute in list_of_metabolite_attributes:
            attributes[attribute] = {}
            for s in network.metabolites:
                if attribute in s.annotation:
                    attributes[attribute][str(s.annotation.get(attribute))] = re.sub(str(s.compartment)+"$","",s.id)
        return attributes

    def attributes_kegg (network, attribute_dictionary): # Specially designed for DB as metabolite name is the KEGG id
        attribute_dictionary['kegg.compound'] = {}
        for s in network.metabolites:
            sth_kegg = re.sub("_"+str(s.compartment)+"$","",s.id)
            sth_id = re.sub(str(s.compartment)+"$","",s.id)          
            attribute_dictionary['kegg.compound'][sth_kegg] = sth_id
        return attribute_dictionary      

    list_of_metabolite_attributes_n1 = attributes (network_1) # List of all the metabolite attributes in network_1
    list_of_metabolite_attributes_n2 = attributes (network_2) # List of all the metabolite attributes in network_2
    list_of_metabolite_attributes_n2 = attributes_kegg (network_2, list_of_metabolite_attributes_n2) # List of all the metabolite attributes in network_2
   

    # Unique metabolites and compartments
    def unique_metabolites_non_unique_metabolites_and_compartments (network):
        unique_metabolites_and_compartments = {}
        for m in network.metabolites:
            mth_unique_metabolite = re.sub(str(m.compartment)+"$","",m.id)
            if not mth_unique_metabolite in unique_metabolites_and_compartments.keys():
                unique_metabolites_and_compartments[mth_unique_metabolite] = {} #{'compartment', 'non_unique_metabolites'}
                unique_metabolites_and_compartments[mth_unique_metabolite]['compartment'] = []
                unique_metabolites_and_compartments[mth_unique_metabolite]['non_unique_metabolites'] = []
            unique_metabolites_and_compartments[mth_unique_metabolite]['compartment'].append(m.compartment)
            unique_metabolites_and_compartments[mth_unique_metabolite]['non_unique_metabolites'].append(m.id)
        return unique_metabolites_and_compartments

    unique_metabolites_and_compartments_n1 = unique_metabolites_non_unique_metabolites_and_compartments(network_1)
    unique_metabolites_and_compartments_n2 = unique_metabolites_non_unique_metabolites_and_compartments(network_2)


    # Equivalent metabolites between n1 and n2
    # 1. by id

    def equivalent_metabolite_by_attribute (list_of_metabolite_attributes_n1, list_of_metabolite_attributes_n2, n1_n2_equivalent_unique_metabolites, attribute):
        n2_list_of_metabolites_with_attribute = list_of_metabolite_attributes_n2[attribute]
        if attribute in list_of_metabolite_attributes_n1.keys():
            for a in n2_list_of_metabolites_with_attribute:
                if a in list_of_metabolite_attributes_n1[attribute].keys():
                    new_equivalence = [list_of_metabolite_attributes_n1[attribute][a],list_of_metabolite_attributes_n2[attribute][a]]
                    if not new_equivalence in n1_n2_equivalent_unique_metabolites:
                        n1_n2_equivalent_unique_metabolites.append(new_equivalence)
        return n1_n2_equivalent_unique_metabolites

    n1_n2_equivalent_unique_metabolites = [] # dictionary created here and used when upgrading reactions, the metabolite ids in the source reactions (from n2) will be changed to fit with the annotation in the target network (n1)
    metabolite_attributes = ['kegg.compound', 'lipidmaps', 'chebi.compound', 'inchi', 'pubchem.compound']
    for attribute in metabolite_attributes:
        n1_n2_equivalent_unique_metabolites = equivalent_metabolite_by_attribute (list_of_metabolite_attributes_n1, list_of_metabolite_attributes_n2, n1_n2_equivalent_unique_metabolites, attribute)
 
    remaining_metabolites_n1 = [x for x in list(unique_metabolites_and_compartments_n1.keys()) if x not in list(set([x[0] for x in n1_n2_equivalent_unique_metabolites]))]
    remaining_metabolites_n2 = [x for x in list(unique_metabolites_and_compartments_n2.keys()) if x not in list(set([x[1] for x in n1_n2_equivalent_unique_metabolites]))]

    # 2. by text similarlity
    for m2 in remaining_metabolites_n2:
        m2th_non_unique_id = unique_metabolites_and_compartments_n2[m2]['non_unique_metabolites'][0]
        for m1 in remaining_metabolites_n1:
            m1th_non_unique_id = unique_metabolites_and_compartments_n1[m1]['non_unique_metabolites'][0]
            # Compare m2th and m1th
            m2th_vs_m1th = compare_metabolite_groups_2 (m1th_non_unique_id,m2th_non_unique_id,network_1,network_2)
            if m2th_vs_m1th == 1:
                print(m1,m2)
                n1_n2_equivalent_unique_metabolites.append(m1,m2)
                remaining_metabolites_n1 = [x for x in remaining_metabolites_n1 if x != m1]
                break

    # Enrich network 1 with the equivalent metabolites in network 2
    for x in n1_n2_equivalent_unique_metabolites: #enrich overlaped metabolites in network_1
        try:
            source_metabolite = network_2.metabolites.get_by_id(unique_metabolites_and_compartments_n2[x[1]]['non_unique_metabolites'][0])
            new_compartments = [c for c in unique_metabolites_and_compartments_n2[x[1]]['compartment'] if not c in unique_metabolites_and_compartments_n1[x[0]]['compartment']]      
            
            for y in unique_metabolites_and_compartments_n1[x[0]]['non_unique_metabolites']: # update the overlaped metabolites
                for a in network_1.metabolites.get_by_id(y).annotation.keys(): # annotation
                    if a in source_metabolite.annotation.keys():
                        network_1.metabolites.get_by_id(y).annotation[a] = source_metabolite.annotation[a]
                #network_1.metabolites.get_by_id(y).charge = copy.deepcopy(source_metabolite.charge)
                #network_1.metabolites.get_by_id(y).elements = copy.deepcopy(source_metabolite.elements)
                #network_1.metabolites.get_by_id(y).formula = copy.deepcopy(source_metabolite.formula)
                network_1.metabolites.get_by_id(y).name = copy.deepcopy(source_metabolite.name)
                network_1.metabolites.get_by_id(y).notes = copy.deepcopy(source_metabolite.notes)

            template_metabolite = copy.deepcopy(network_1.metabolites.get_by_id(y))

            for z in new_compartments: # add metabolites in compartments not annotated in the target network 
                zth_metabolite = copy.deepcopy(template_metabolite)
                zth_metabolite.id = re.sub(zth_metabolite.compartment+'$', z, zth_metabolite.id)
                zth_metabolite.compartment = z 
                network_1.add_metabolites(zth_metabolite)
        except:
            continue

    # Add the non-overlapped metabolites
    # remaining_metabolites_n2 = [x for x in list(unique_metabolites_and_compartments_n2.keys()) if x not in list(set([x[1] for x in n1_n2_equivalent_unique_metabolites]))]

    for x in remaining_metabolites_n2:
        non_unique_metabolites = unique_metabolites_and_compartments_n2[x]['non_unique_metabolites']
        for y in non_unique_metabolites:
            xyth_metabolite = copy.deepcopy(network_2.metabolites.get_by_id(y))
            network_1.add_metabolites(xyth_metabolite)

    return network_1, n1_n2_equivalent_unique_metabolites


# Step 2: transter genes from source network to target network
'Merge gprs from two different reactions by logic "or"'
def merge_gprs_old(reaction1, reaction2):
    """
    For two reactions, the function takes the GPRs removes duplicate items (enzymes, isoenzymes and complexes) and join the two GPRs using logic "or"

    Inputs
    ----------
    reaction1 and reaction2 : metabolic reaction in Cobra format including gpr attribute
    
    Output
    -------
    joint_gpr: join GPR as cobra.core.gene.GPR object 
    """
    gpr1 = [x.strip() for x in reaction1.gpr.to_string().split("or")]
    gpr2 = [x.strip() for x in reaction2.gpr.to_string().split("or")]
    join_gpr = cobra.core.gene.GPR.from_string(re.sub('@', ' ', re.sub(' ', ' or ', re.sub(' and ', '@and@',' '.join(list(set(gpr1 + gpr2))).strip()))))
    return join_gpr


def divide_gpr_in_ors(ast) -> List:
    """Divide GPR in AST form recursively until a Gene or a AND operation is found."""
    if isinstance(ast, Name) or isinstance(ast, str):
        # base case 1: a gene
        return [ast]
    elif isinstance(ast, BoolOp):
        if isinstance(ast.op, Or):
            return reduce(lambda x, y: x + divide_gpr_in_ors(y), ast.values, [])
        if isinstance(ast.op, And):
            # base case 2: AND is returned as is
            return [ast]


def compare_ast(ast_left, ast_right) -> bool:
    """Compare recursively that 2 AST of genes are the same."""
    # base case 1: the asts are str or Name -> str comparison
    if isinstance(ast_left, str) and isinstance(ast_right, str):
        if ast_left != ast_right:
            return False
    elif isinstance(ast_left, Name) and isinstance(ast_right, Name):
        if ast_left.id != ast_right.id:
            return False
    elif isinstance(ast_left, BoolOp) and isinstance(ast_right, BoolOp):
        for left_leaf in ast_left.values:
            if not any(
                compare_ast(left_leaf, leaf_right) for leaf_right in ast_right.values
            ):
                return False
    else:
        # base case 2: types do not match
        return False
    return True


def merge_ast_gprs(ast_left: Expression, ast_right: Expression) -> GPR:
    """Merge two GPRs represented as GPR `ast.Expression`."""
    or_left = divide_gpr_in_ors(ast_left)
    or_right = divide_gpr_in_ors(ast_right)
    for left_node in or_left:
        if any(compare_ast(left_node, right_node) for right_node in or_right):
            pass
        else:
            or_right.append(left_node)
    return cobra.core.gene.GPR(Expression(BoolOp(Or(), or_right)))


def merge_gprs(
    reaction_left: Reaction, reaction_right: Reaction
) -> GPR:
    """Merge the GPRs of two reactions at the AST level.

    Algorithm:

    1. ast_left, ast_right <- Parse AST (OR and ANDs) for reaction 1 and reaction 2.
    2. ors_left, ors_right <- Divide ASTs in Ors (the resulting nodes are AND or Genes).
    3. for left_node in ors_left:
        if left_node not found in ors_right:
            ors_right += left_node
    4. return GPR(ors_right)

    """
    # this is divided in merge_ast_gprs for easier testing
    return merge_ast_gprs(
        ast_parse(reaction_left.gene_reaction_rule).body[0].value,
        ast_parse(reaction_right.gene_reaction_rule).body[0].value,
    )


'Enrich gene list and annotation in network_1 from network_2'
def network_genes_merge (network_1, network_2):
    """
    Genes from two metabolic models are compared based on ID. If a match is found the information from network_2 replaces the infromation in network_1

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with new metabolites from network_2
    """
    [network_1.genes.add(x) for x in network_2.genes if x not in network_1.genes]
    return network_1


def network_genes_merge_2 (network_1, network_2):
    """
    Genes from two metabolic models are compared based on ID. If a match is found the information from network_2 replaces the infromation in network_1

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with new metabolites from network_2
    n1_n2_equivalent_genes : list of equivalent genes between network_1 and network_2
    """
    # Attributes
    def attributes_genes (network):
        list_of_gene_attributes = [] # List of all the metabolite attributes in network_1
        for s in network.genes:
            list_of_gene_attributes = list_of_gene_attributes + [x for x in s.annotation]
        list_of_gene_attributes = list(set(list_of_gene_attributes)) 
        attributes = {} # dictionary containing all the attributes of all the metabolites in network_1
        for attribute in list_of_gene_attributes:
            if attribute != 'sbo' and attribute != 'hgcn.symbol' and attribute != 'hgnc.symbol':
                attributes[attribute] = {}
                for g in network.genes:
                    if attribute in g.annotation.keys():
                        if isinstance(g.annotation[attribute], list):
                            for ss in g.annotation[attribute]:
                                attributes[attribute][ss] = g.id
                        else:
                            attributes[attribute][g.annotation[attribute]] = g.id
        return attributes
    
    list_of_gene_attributes_n1 = attributes_genes (network_1) # List of all the metabolite attributes in network_1
    list_of_gene_attributes_n2 = attributes_genes (network_2) # List of all the metabolite attributes in network_2 

    genes_attributes = list(set(list(list_of_gene_attributes_n1.keys())+list(list_of_gene_attributes_n2.keys())))


    # Equivalent genes between n1 and n2
    # 1. by id
    def equivalent_genes_by_attribute (list_of_gene_attributes_n1, list_of_gene_attributes_n2, n1_n2_equivalent_genes, attribute):
        n2_list_of_genes_with_attribute = list_of_gene_attributes_n2[attribute]
        if attribute in list_of_gene_attributes_n1.keys():
            for a in n2_list_of_genes_with_attribute:
                if a in list_of_gene_attributes_n1[attribute].keys():
                    new_equivalence = [list_of_gene_attributes_n1[attribute][a],list_of_gene_attributes_n2[attribute][a]]
                    if not new_equivalence in n1_n2_equivalent_genes:
                        n1_n2_equivalent_genes.append(new_equivalence)
        return n1_n2_equivalent_genes

    n1_n2_equivalent_genes = [] # dictionary created here and used when upgrading reactions, the metabolite ids in the source reactions (from n2) will be changed to fit with the annotation in the target network (n1)
    for attribute in genes_attributes:
        n1_n2_equivalent_genes = equivalent_genes_by_attribute (list_of_gene_attributes_n1, list_of_gene_attributes_n2, n1_n2_equivalent_genes, attribute)
        
    remaining_genes_n2 = [x.id for x in network_2.genes if x.id not in list(set([x[1] for x in n1_n2_equivalent_genes]))]

    # Enrich network 1 with the equivalent genes in network 2
    for x in n1_n2_equivalent_genes: #enrich overlaped genes in network_1
        source_gene = network_2.genes.get_by_id(x[1])
        for a in network_1.genes.get_by_id(x[0]).annotation.keys(): # annotation
            if a != 'sbo'and a != 'hgcn.symbol' and a != 'hgnc.symbol' and a in source_gene.annotation.keys():
                network_1.genes.get_by_id(x[0]).annotation[a] = source_gene.annotation[a]

   # Add the non-overlapped metabolites
    for x in remaining_genes_n2:
        xth_gene = copy.deepcopy(network_2.genes.get_by_id(x))
        network_1.genes.add(xth_gene)

    return network_1, n1_n2_equivalent_genes


# Step 3: Merge and add reactions from network 2 to network 1
'Enrich reaction list and annotation in network_1 from network_2'
def network_reactions_merge (network_1, network_2, n1_n2_equivalent_metabolites, equivalent_genes):
    """
    Reactions from two metabolic models are compared based on the associated metabolites. If a match is found the information from network_2 replaces the infromation in network_1 or new reactions from network_2 are added to network_1 if no match is found

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    n1_n2_equivalent_metabolites : dictionary with equivalent metabolite ID between network_1 and network_2
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated reaction information from network_2 for the overlaped reactions and new reactions from network_2 in case they were not annotated in network_1 previously
    """
    reaction_equivalent = []
    reaction_no_equivalent = []

    h_in_n_1 = [x.id for x in network_1.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_1
    h_in_n_2 = [x.id for x in network_2.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_2
    reactions_overlaped_n1_n2 = [] # list of checked reactions in n_1
    reaction_added_from_n2_to_n1 = [] # list of reactions in n_2 added to n_1
    # equivalent_metabolites_key_list = [x[1] for x in n1_n2_equivalent_metabolites]
    equivalent_metabolites_key_list = {x[1]:x[0] for x in n1_n2_equivalent_metabolites} # key:n2 id, value: n1 id
    n2_n1_equivalent_genes = {x[1]:x[0] for x in equivalent_genes} # key:n2 id, value: n1 id

    #equivalent_metabolites_key_list = [x for x in n1_n2_equivalent_metabolites.keys()]
    for reaction_n2 in network_2.reactions:
        # 2.1: Identify metabolites in reaction from network 2 
        metabolites_in_n2 = [x.id for x in reaction_n2.metabolites]
        reaction = reaction_n2.reaction
        for x in metabolites_in_n2:
            x_unique = re.sub(network_2.metabolites.get_by_id(x).compartment+'$', '', network_2.metabolites.get_by_id(x).id)
            if x_unique in equivalent_metabolites_key_list.keys():
                equivalent_metabolites_key_list.keys()      
                reaction = re.sub(x, equivalent_metabolites_key_list[x_unique]+network_2.metabolites.get_by_id(x).compartment,reaction_n2.reaction)
                reaction_n2.reaction = reaction    
                reaction_n2.metabolites.update()
        if len(reaction_n2.reactants) > 0:
            substrates_n2 = [x.id for x in reaction_n2.reactants if x.id not in h_in_n_2]
        if len(reaction_n2.products) > 0:
            products_n2 = [x.id for x in reaction_n2.products if x.id not in h_in_n_2]        
        for reaction_n1 in network_1.reactions:
            # 2.2: Identify metabolites in reaction from network 1
            if reaction_n1.id not in reactions_overlaped_n1_n2 and reaction_n2.compartments == reaction_n1.compartments: # only compare if reaction_n1 was not matched before and both reactions are in the same compartments
                if len(reaction_n1.reactants) > 0: 
                    substrates_n1 = [x.id for x in reaction_n1.reactants if x.id not in h_in_n_1]
                if len(reaction_n1.products) > 0: 
                    products_n1 = [x.id for x in reaction_n1.products if x.id not in h_in_n_1]
                # 2.3: Jaccard Index comparison of substrates and products of reactions from network 1 and 2
                if len(reaction_n2.products) == 0 and len(reaction_n1.products) == 0: # reaction in network 1 and 2 are exchange reactions
                    # Jaccard Index
                    test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1))
                    test_redundancy_1[0] = test_redundancy_1[0]*2
                    test_redundancy_2 = list((0,[0],[0]))
                if len(reaction_n2.products) != 0 and len(reaction_n1.products) != 0: # reaction in network 1 and 2 are not exchange reactions
                    # Jaccard Index
                    test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1) + compare_metabolite_groups(products_n2, products_n1, network_2, network_1))
                    test_redundancy_1[0] = test_redundancy_1[0]+test_redundancy_1[2]
                    del test_redundancy_1[2]
                    reversibility_factor = reaction_n2.reversibility*reaction_n1.reversibility # 1: both reactions are reversible, 0: one or both reactions are irreversible
                    test_redundancy_2 = list(compare_metabolite_groups(substrates_n2, products_n1, network_2, network_1) + compare_metabolite_groups(products_n2, substrates_n1, network_2, network_1))
                    test_redundancy_2[0] = reversibility_factor*(test_redundancy_2[0]+test_redundancy_2[2])
                    del test_redundancy_2[2]               
                # 2.4: Modify/add reaction in network 1 based on the results from 2.3    
                if test_redundancy_1[0] == 2 or test_redundancy_2[0] == 2: # if the reactions are the same replace the information in reaction_n1 by the information in reaction_n2
                    # 2.4.1: Add reaction and metabolite attributes from netwrok 2 to network 1
                    if test_redundancy_1[0] < test_redundancy_2[0]: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible               
                    # if sum(test_redundancy_1[1]+test_redundancy_1[2]) < sum(test_redundancy_2[1]+test_redundancy_2[2])*reversibility_factor: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible
                        reaction = reaction.split(' <=> ')[1]+' <=> '+reaction.split(' <=> ')[0]
                    reaction_n1.reaction = copy.deepcopy(reaction)
                    reaction_n1.metabolites.update()
                    # 2.4.2: Add gene attributes from netwrok 2 to network 1              
                    # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
                    gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
                    gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
                    genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
                    for g in genes_in_n2:
                        if g in n2_n1_equivalent_genes.keys():                 
                            gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                            new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                            reaction_n1.gpr = new_gpr
                            reaction_n1.update_genes_from_gpr()
                    # 2.4.2: Add other attributes from netwrok 2 to network 1
                    reaction_n1.annotation = copy.deepcopy(reaction_n2.annotation)
                    reaction_n1.name = copy.deepcopy(reaction_n2.name)
                    reactions_overlaped_n1_n2.append(reaction_n1.id) # Add the reaction_n1 to the list of matched reaction in n1
                    reaction_equivalent = reaction_equivalent.append([reaction_n1.id,reaction_n2.id])
                    reaction_added_from_n2_to_n1.append(reaction_n2)
                    break # once matched the inner loop can stop 
    # Add reactions from n2 not matching with any reaction in n1    

    for x in network_2.reactions:
        if x not in reaction_added_from_n2_to_n1:
            print(x)
            # Give a new id consistent with network1 annotation
            if x.id in network_1.reactions: # in case the id, btt chance, already exist
                last_reaction = sorted([y.id for y in network_1.reactions],reverse=True)[0]
                next_reaction = re.findall("^([A-Z,a-z]*?)[0-9]",last_reaction)[0]+str(int((re.findall("^[A-Z,a-z]*?([0-9].*?)$",last_reaction)[0]))+1)
                x.id = next_reaction

            # Change metabolites id according to network_1
            metabolites_in_n2 = [y.id for y in x.metabolites]
            reaction = x.reaction
            for y in metabolites_in_n2:
                if not network_2.metabolites.get_by_id(y).compartment:
                    metabolite_compartment = network_1.metabolites.get_by_id(y).compartment
                else:
                    metabolite_compartment = network_2.metabolites.get_by_id(y).compartment
                y_unique = re.sub(metabolite_compartment+'$', '', network_2.metabolites.get_by_id(y).id)
                if y_unique in equivalent_metabolites_key_list.keys():  
                    reaction = re.sub(y, equivalent_metabolites_key_list[y_unique]+network_2.metabolites.get_by_id(y).compartment,x.reaction)
                    x.reaction = reaction    
                    x.metabolites.update()

            # Change gene id according to network_1

            gene_rules = copy.deepcopy(x.gene_reaction_rule)
            gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
            genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
            for g in genes_in_n2:
                if g in n2_n1_equivalent_genes.keys():                 
                    gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                    new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                    x.gpr = new_gpr
                    x.update_genes_from_gpr()

            # Add the reaction to network_1

            network_1.add_reactions([x])

            reaction_no_equivalent.append([x.id])
    return network_1, reaction_equivalent, reaction_no_equivalent


def network_reactions_merge_2 (network_1, network_2, n1_n2_equivalent_metabolites, equivalent_genes): # deprecated until I find out why the models cannot generate biomass
    """
    Reactions from two metabolic models are compared based on the associated metabolites. If a match is found the information from network_2 replaces the infromation in network_1 or new reactions from network_2 are added to network_1 if no match is found

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    n1_n2_equivalent_metabolites : dictionary with equivalent metabolite ID between network_1 and network_2
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated reaction information from network_2 for the overlaped reactions and new reactions from network_2 in case they were not annotated in network_1 previously
    overlap_reactions: list
    not_overlap_reactions: list
    """

    ''' 1. Warm up'''
    from collections import Counter

    tolerance = 1.75 # values between 0 and 2
    h_in_n_1 = [x.id for x in network_1.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_1
    h_in_n_2 = [x.id for x in network_2.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_2

    # Evaluate reactions in network_2
    internal_rxns = con_helpers.get_internals(network_2)
    imbalanced_reactions = get_ids(consistency.find_mass_unbalanced_reactions(internal_rxns))
    ListOfMetFrom = {x.id: x.formula for x in network_2.metabolites}
    imbalanced_reactions_network_2 = [(network_2.reactions.get_by_id(x).id, imbalance_test (network_2.reactions.get_by_id(x), ListOfMetFrom)) for x in imbalanced_reactions]


    ''' 2. Identify Equivalent Reactions '''

    print('Matching reactions')

    # 2.1. By E.C.
    ec_network_2 = [(x.id, x.annotation['ec-code']) for x in network_2.reactions if 'ec-code' in x.annotation.keys()]
    ec_network_1 = [(x.id, x.annotation['ec-code']) for x in network_1.reactions if 'ec-code' in x.annotation.keys()]
    ec_overlap = list(sorted([[x[0],y[0]] for x in ec_network_2 for y in ec_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.2. By KEGG
    kegg_network_2 = [(x.id, x.annotation['kegg.reaction']) for x in network_2.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_network_1 = [(x.id, x.annotation['kegg.reaction']) for x in network_1.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_overlap = list(sorted([[x[0],y[0]] for x in kegg_network_2 for y in kegg_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.3. Combine KEGG + E.C.
    from more_itertools import locate

    def find_indices(list_to_check, item_to_find):
        indices = locate(list_to_check, lambda x: x == item_to_find)
        return list(indices)

    kegg_ec_overlap = [x for x in kegg_overlap if x in ec_overlap]
    counts = dict(Counter([x[0] for x in kegg_ec_overlap]))
    dm = {key:value for key, value in counts.items() if value > 1}
    uniques_matches = [x for x in kegg_ec_overlap if not x[0] in dm.keys()]
    duplicate_matches = [(x,[kegg_ec_overlap[y][1] for y in find_indices([z[0]for z in kegg_ec_overlap], x)]) for x in dm.keys()]

    # 2.4. Find equivalent reactions between networks

    perfect_matches = []

    for x in uniques_matches+duplicate_matches: # When more than one equivalence from ec+kegg, choose the best

        xth_n2 = x[0]
        xth_n1 = x[1]
        if isinstance(xth_n1,str): xth_n1 = [xth_n1]

        reaction_n2 = network_2.reactions.get_by_id(xth_n2)
        if len(reaction_n2.reactants) > 0:
            substrates_n2 = [x.id for x in reaction_n2.reactants if x.id not in h_in_n_2]
        if len(reaction_n2.products) > 0:
            products_n2 = [x.id for x in reaction_n2.products if x.id not in h_in_n_2]  

        test_redundancy_ref = 0
        best_match = []

        for y in xth_n1:
            reaction_n1 = network_1.reactions.get_by_id(y)
            if len(reaction_n1.reactants) > 0: 
                substrates_n1 = [x.id for x in reaction_n1.reactants if x.id not in h_in_n_1]
            if len(reaction_n1.products) > 0: 
                products_n1 = [x.id for x in reaction_n1.products if x.id not in h_in_n_1]

            if len(reaction_n2.products) == 0 and len(reaction_n1.products) == 0: # reaction in network 1 and 2 are exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]*2
                test_redundancy_2 = list((0,[0],[0]))
            if len(reaction_n2.products) != 0 and len(reaction_n1.products) != 0: # reaction in network 1 and 2 are not exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1) + compare_metabolite_groups(products_n2, products_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]+test_redundancy_1[2]
                del test_redundancy_1[2]
                reversibility_factor = reaction_n2.reversibility*reaction_n1.reversibility # 1: both reactions are reversible, 0: one or both reactions are irreversible
                test_redundancy_2 = list(compare_metabolite_groups(substrates_n2, products_n1, network_2, network_1) + compare_metabolite_groups(products_n2, substrates_n1, network_2, network_1))
                test_redundancy_2[0] = reversibility_factor*(test_redundancy_2[0]+test_redundancy_2[2])
                del test_redundancy_2[2]    
            if max(test_redundancy_1[0],test_redundancy_2[0]) >= tolerance: # if the reactions are the same replace add the metabolite to the list of equivalent
                if max(test_redundancy_1[0],test_redundancy_2[0]) >= test_redundancy_ref:
                    test_redundancy_ref = max(test_redundancy_1[0],test_redundancy_2[0])
                    best_match = [xth_n2, y,test_redundancy_1[0], test_redundancy_2[0]]
            if best_match:
                perfect_matches.append(best_match)
                best_match = []


    ''' 3. Change ids and annotation in network 2 to be compatible with network 1 '''

    print('Making annotation compatible between networks')

    equivalent_metabolites_key_list = {x[1]:x[0] for x in n1_n2_equivalent_metabolites} # key:n2 id, value: n1 id
    n2_n1_equivalent_genes = {x[1]:x[0] for x in equivalent_genes} # key:n2 id, value: n1 id

    for reaction_n2 in network_2.reactions:
        # 3.1: Identify metabolites in reaction from network 2 
        metabolites_in_n2 = [x.id for x in reaction_n2.metabolites]
        reaction = reaction_n2.reaction
        for x in metabolites_in_n2:
            x_unique = re.sub(network_2.metabolites.get_by_id(x).compartment+'$', '', network_2.metabolites.get_by_id(x).id)
            if x_unique in equivalent_metabolites_key_list.keys():
                #equivalent_metabolites_key_list.keys()      
                reaction = re.sub(x, equivalent_metabolites_key_list[x_unique]+network_2.metabolites.get_by_id(x).compartment,reaction_n2.reaction)
                reaction_n2.reaction = reaction    
                reaction_n2.metabolites.update()
        # 3.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                reaction_n2.gpr = new_gpr
                reaction_n2.update_genes_from_gpr()


    ''' 4. Enrich network 1 with the overlaped reactions '''

    # 4.1: Add reaction and metabolite attributes from overlap reaction in netwrok 2 to network 1
    print('Merging overlap reactions')
    for x in perfect_matches:
        reaction_n1 = network_1.reactions.get_by_id(x[1])
        reaction_n2 = network_2.reactions.get_by_id(x[0])

        # 2.4.1: Add reaction and metabolite attributes from netwrok 2 to network 1
        reaction = reaction_n2.reaction
        if x[2] < x[3]: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible               
            reaction = reaction.split(' <=> ')[1]+' <=> '+reaction.split(' <=> ')[0]
        if reaction_n1.reversibility == True:
            rev_n2 = re.findall('(...)>', reaction_n2.reaction)[0].strip()+'>'
            reaction = re.sub(rev_n2, '<=>', reaction) #Impose reversibility from network 1

        reaction_n1.reaction = copy.deepcopy(reaction)
        reaction_n1.metabolites.update()
        
        # 2.4.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))

        reaction_n1.gpr = copy.deepcopy(new_gpr)
        reaction_n1.update_genes_from_gpr()
                        
        # 2.4.2: Add other attributes from netwrok 2 to network 1
        # annotation
        new_annotation = copy.deepcopy(reaction_n1.annotation)
        annotation_overlap = list(reduce(lambda i, j: i & j, (set(x.annotation.keys()) for x in [reaction_n1]+[reaction_n2]))) # common annotations between reaction_n1 and reaction_n2
        for a in annotation_overlap:
            new_annotation[a] = reaction_n2.annotation[a] 
        for a_n2 in reaction_n2.annotation.keys():
            new_annotation[a_n2] = reaction_n2.annotation[a_n2] 

        reaction_n1.annotation = copy.deepcopy(new_annotation)

        #name
        reaction_n1.name = copy.deepcopy(reaction_n2.name)


    # 4.2: Add reaction and metabolite attributes from not overlap reaction in netwrok 2 to network 1
    print('Adding new reactions')
    reactions_n2_not_overlap = [x for x in network_2.reactions if not x.id in [y[0] for y in perfect_matches] and not x.id in imbalanced_reactions]
    
    for x in reactions_n2_not_overlap:
        network_1.reactions.add(x)
        network_1.reactions.get_by_id(x.id).metabolites.update()


    ''' 5. Summary '''
    overlap_reactions = [y[0] for y in perfect_matches]
    not_overlap_reactions = [x.id for x in network_2.reactions if not x.id in [y[0] for y in perfect_matches]]

    return network_1, overlap_reactions, not_overlap_reactions


def network_reactions_merge_3 (network_1, network_2, n1_n2_equivalent_metabolites, equivalent_genes): # like network_reactions_merge but it is added a consistency test before adding each reaction
    """
    Reactions from two metabolic models are compared based on the associated metabolites. If a match is found the information from network_2 replaces the infromation in network_1 or new reactions from network_2 are added to network_1 if no match is found

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    n1_n2_equivalent_metabolites : dictionary with equivalent metabolite ID between network_1 and network_2
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated reaction information from network_2 for the overlaped reactions and new reactions from network_2 in case they were not annotated in network_1 previously
    """
    reaction_equivalent = []
    reaction_no_equivalent = []

    h_in_n_1 = [x.id for x in network_1.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_1
    h_in_n_2 = [x.id for x in network_2.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_2
    reactions_overlaped_n1_n2 = [] # list of checked reactions in n_1
    reaction_added_from_n2_to_n1 = [] # list of reactions in n_2 added to n_1
    # equivalent_metabolites_key_list = [x[1] for x in n1_n2_equivalent_metabolites]
    equivalent_metabolites_key_list = {x[1]:x[0] for x in n1_n2_equivalent_metabolites} # key:n2 id, value: n1 id
    n2_n1_equivalent_genes = {x[1]:x[0] for x in equivalent_genes} # key:n2 id, value: n1 id

    #equivalent_metabolites_key_list = [x for x in n1_n2_equivalent_metabolites.keys()]
    for reaction_n2 in network_2.reactions:
        # 2.1: Identify metabolites in reaction from network 2 
        metabolites_in_n2 = [x.id for x in reaction_n2.metabolites]
        reaction = reaction_n2.reaction
        for x in metabolites_in_n2:
            x_unique = re.sub(network_2.metabolites.get_by_id(x).compartment+'$', '', network_2.metabolites.get_by_id(x).id)
            if x_unique in equivalent_metabolites_key_list.keys():
                equivalent_metabolites_key_list.keys()      
                reaction = re.sub(x, equivalent_metabolites_key_list[x_unique]+network_2.metabolites.get_by_id(x).compartment,reaction_n2.reaction)
                reaction_n2.reaction = reaction    
                reaction_n2.metabolites.update()
        if len(reaction_n2.reactants) > 0:
            substrates_n2 = [x.id for x in reaction_n2.reactants if x.id not in h_in_n_2]
        if len(reaction_n2.products) > 0:
            products_n2 = [x.id for x in reaction_n2.products if x.id not in h_in_n_2]        
        for reaction_n1 in network_1.reactions:
            # 2.2: Identify metabolites in reaction from network 1
            if reaction_n1.id not in reactions_overlaped_n1_n2 and reaction_n2.compartments == reaction_n1.compartments: # only compare if reaction_n1 was not matched before and both reactions are in the same compartments
                if len(reaction_n1.reactants) > 0: 
                    substrates_n1 = [x.id for x in reaction_n1.reactants if x.id not in h_in_n_1]
                if len(reaction_n1.products) > 0: 
                    products_n1 = [x.id for x in reaction_n1.products if x.id not in h_in_n_1]
                # 2.3: Jaccard Index comparison of substrates and products of reactions from network 1 and 2
                if len(reaction_n2.products) == 0 and len(reaction_n1.products) == 0: # reaction in network 1 and 2 are exchange reactions
                    # Jaccard Index
                    test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1))
                    test_redundancy_1[0] = test_redundancy_1[0]*2
                    test_redundancy_2 = list((0,[0],[0]))
                if len(reaction_n2.products) != 0 and len(reaction_n1.products) != 0: # reaction in network 1 and 2 are not exchange reactions
                    # Jaccard Index
                    test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1) + compare_metabolite_groups(products_n2, products_n1, network_2, network_1))
                    test_redundancy_1[0] = test_redundancy_1[0]+test_redundancy_1[2]
                    del test_redundancy_1[2]
                    reversibility_factor = reaction_n2.reversibility*reaction_n1.reversibility # 1: both reactions are reversible, 0: one or both reactions are irreversible
                    test_redundancy_2 = list(compare_metabolite_groups(substrates_n2, products_n1, network_2, network_1) + compare_metabolite_groups(products_n2, substrates_n1, network_2, network_1))
                    test_redundancy_2[0] = reversibility_factor*(test_redundancy_2[0]+test_redundancy_2[2])
                    del test_redundancy_2[2]               
                # 2.4: Modify/add reaction in network 1 based on the results from 2.3    
                if test_redundancy_1[0] == 2 or test_redundancy_2[0] == 2: # if the reactions are the same replace the information in reaction_n1 by the information in reaction_n2
                    # 2.4.1: Add reaction and metabolite attributes from netwrok 2 to network 1
                    if test_redundancy_1[0] < test_redundancy_2[0]: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible               
                    # if sum(test_redundancy_1[1]+test_redundancy_1[2]) < sum(test_redundancy_2[1]+test_redundancy_2[2])*reversibility_factor: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible
                        reaction = reaction.split(' <=> ')[1]+' <=> '+reaction.split(' <=> ')[0]
                    reaction_n1.reaction = copy.deepcopy(reaction)
                    reaction_n1.metabolites.update()
                    # 2.4.2: Add gene attributes from netwrok 2 to network 1              
                    # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
                    gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
                    gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
                    genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
                    for g in genes_in_n2:
                        if g in n2_n1_equivalent_genes.keys():                 
                            gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                            new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                            reaction_n1.gpr = new_gpr
                            reaction_n1.update_genes_from_gpr()
                    # 2.4.2: Add other attributes from netwrok 2 to network 1
                    reaction_n1.annotation = copy.deepcopy(reaction_n2.annotation)
                    reaction_n1.name = copy.deepcopy(reaction_n2.name)
                    reactions_overlaped_n1_n2.append(reaction_n1.id) # Add the reaction_n1 to the list of matched reaction in n1
                    reaction_equivalent = reaction_equivalent.append([reaction_n1.id,reaction_n2.id])
                    reaction_added_from_n2_to_n1.append(reaction_n2)
                    break # once matched the inner loop can stop 
    # Add reactions from n2 not matching with any reaction in n1    

    for x in network_2.reactions:
        if x not in reaction_added_from_n2_to_n1:
            network_1_test = copy.deepcopy(network_1)
            print(x)
            # Give a new id consistent with network1 annotation
            if x.id in network_1_test.reactions: # in case the id, btt chance, already exist
                last_reaction = sorted([y.id for y in network_1_test.reactions],reverse=True)[0]
                next_reaction = re.findall("^([A-Z,a-z]*?)[0-9]",last_reaction)[0]+str(int((re.findall("^[A-Z,a-z]*?([0-9].*?)$",last_reaction)[0]))+1)
                x.id = next_reaction

            # Change metabolites id according to network_1_test
            metabolites_in_n2 = [y.id for y in x.metabolites]
            reaction = x.reaction
            for y in metabolites_in_n2:
                if not network_2.metabolites.get_by_id(y).compartment:
                    metabolite_compartment = network_1_test.metabolites.get_by_id(y).compartment
                else:
                    metabolite_compartment = network_2.metabolites.get_by_id(y).compartment
                y_unique = re.sub(metabolite_compartment+'$', '', network_2.metabolites.get_by_id(y).id)
                if y_unique in equivalent_metabolites_key_list.keys():  
                    reaction = re.sub(y, equivalent_metabolites_key_list[y_unique]+network_2.metabolites.get_by_id(y).compartment,x.reaction)
                    x.reaction = reaction    
                    x.metabolites.update()

            # Change gene id according to network_1_test

            gene_rules = copy.deepcopy(x.gene_reaction_rule)
            gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
            genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
            for g in genes_in_n2:
                if g in n2_n1_equivalent_genes.keys():                 
                    gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                    new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                    x.gpr = new_gpr
                    x.update_genes_from_gpr()

            # Add the reaction to network_1_test

            network_1_test.add_reactions([x])

            model_consistency = consistency.check_stoichiometric_consistency(network_1_test)
            
            if model_consistency == True: network_1 = copy.deepcopy(network_1_test)

            reaction_no_equivalent.append([x.id])
    return network_1, reaction_equivalent, reaction_no_equivalent


def network_reactions_merge_4 (network_1, network_2, n1_n2_equivalent_metabolites, equivalent_genes): # combines version 2 and 4: fast overlap search + robust addition of new reaction with a testing step
    """
    Reactions from two metabolic models are compared based on the associated metabolites. If a match is found the information from network_2 replaces the infromation in network_1 or new reactions from network_2 are added to network_1 if no match is found

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    n1_n2_equivalent_metabolites : dictionary with equivalent metabolite ID between network_1 and network_2
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated reaction information from network_2 for the overlaped reactions and new reactions from network_2 in case they were not annotated in network_1 previously
    overlap_reactions: list
    not_overlap_reactions: list
    """

    ''' 1. Warm up'''
    from collections import Counter

    tolerance = 1.75 # values between 0 and 2
    h_in_n_1 = [x.id for x in network_1.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_1
    h_in_n_2 = [x.id for x in network_2.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_2
    ListOfMetFrom = {x.id: x.formula for x in network_1.metabolites}
    ListOfMetFrom_inverse = {x.formula: x.id for x in network_1.metabolites}
      
    Atom_ID = {'H2O1': 'MAM02040', 'H1': 'MAM02039', 'H2O': 'MAM02040', 'H': 'MAM02039', 'Fe': 'MAM01821','F': 'MAM01821', 'X': 'MAM01823', 'R': 'MAM03573', 'Na': 'MAM02519', 'K': 'MAM02200', 'Ca': 'MAM01413'}
    Atom_ID_Inverse = {'MAM02040': 'H2O1', 'MAM02039': 'H1', 'MAM02040': 'H2O', 'MAM02039': 'H', 'MAM01821': 'Fe','MAM01821': 'F', 'MAM01823': 'X', 'MAM03573': 'R', 'MAM02519': 'Na', 'MAM02200': 'K', 'MAM01413': 'Ca'}

    network_atom_equivalent = {'C00000': 'MAM03573'}

    # Evaluate reactions in network_2
    internal_rxns = con_helpers.get_internals(network_2)
    imbalanced_reactions = get_ids(consistency.find_mass_unbalanced_reactions(internal_rxns))
    list_of_innconsistent_reactions = []

    ''' 2. Identify Equivalent Reactions '''

    print('Matching reactions')

    # 2.1. By E.C.
    ec_network_2 = [(x.id, x.annotation['ec-code']) for x in network_2.reactions if 'ec-code' in x.annotation.keys()]
    ec_network_1 = [(x.id, x.annotation['ec-code']) for x in network_1.reactions if 'ec-code' in x.annotation.keys()]
    ec_overlap = list(sorted([[x[0],y[0]] for x in ec_network_2 for y in ec_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.2. By KEGG
    kegg_network_2 = [(x.id, x.annotation['kegg.reaction']) for x in network_2.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_network_1 = [(x.id, x.annotation['kegg.reaction']) for x in network_1.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_overlap = list(sorted([[x[0],y[0]] for x in kegg_network_2 for y in kegg_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.3. Combine KEGG + E.C.
    from more_itertools import locate

    def find_indices(list_to_check, item_to_find):
        indices = locate(list_to_check, lambda x: x == item_to_find)
        return list(indices)

    kegg_ec_overlap = [x for x in kegg_overlap if x in ec_overlap]
    counts = dict(Counter([x[0] for x in kegg_ec_overlap]))
    dm = {key:value for key, value in counts.items() if value > 1}
    uniques_matches = [x for x in kegg_ec_overlap if not x[0] in dm.keys()]
    duplicate_matches = [(x,[kegg_ec_overlap[y][1] for y in find_indices([z[0]for z in kegg_ec_overlap], x)]) for x in dm.keys()]

    # 2.4. Find equivalent reactions between networks

    perfect_matches = []

    for x in uniques_matches+duplicate_matches: # When more than one equivalence from ec+kegg, choose the best

        xth_n2 = x[0]
        xth_n1 = x[1]
        if isinstance(xth_n1,str): xth_n1 = [xth_n1]

        reaction_n2 = network_2.reactions.get_by_id(xth_n2)
        if len(reaction_n2.reactants) > 0:
            substrates_n2 = [x.id for x in reaction_n2.reactants if x.id not in h_in_n_2]
        if len(reaction_n2.products) > 0:
            products_n2 = [x.id for x in reaction_n2.products if x.id not in h_in_n_2]  

        test_redundancy_ref = 0
        best_match = []

        for y in xth_n1:
            reaction_n1 = network_1.reactions.get_by_id(y)
            if len(reaction_n1.reactants) > 0: 
                substrates_n1 = [x.id for x in reaction_n1.reactants if x.id not in h_in_n_1]
            if len(reaction_n1.products) > 0: 
                products_n1 = [x.id for x in reaction_n1.products if x.id not in h_in_n_1]

            if len(reaction_n2.products) == 0 and len(reaction_n1.products) == 0: # reaction in network 1 and 2 are exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]*2
                test_redundancy_2 = list((0,[0],[0]))
            if len(reaction_n2.products) != 0 and len(reaction_n1.products) != 0: # reaction in network 1 and 2 are not exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1) + compare_metabolite_groups(products_n2, products_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]+test_redundancy_1[2]
                del test_redundancy_1[2]
                reversibility_factor = reaction_n2.reversibility*reaction_n1.reversibility # 1: both reactions are reversible, 0: one or both reactions are irreversible
                test_redundancy_2 = list(compare_metabolite_groups(substrates_n2, products_n1, network_2, network_1) + compare_metabolite_groups(products_n2, substrates_n1, network_2, network_1))
                test_redundancy_2[0] = reversibility_factor*(test_redundancy_2[0]+test_redundancy_2[2])
                del test_redundancy_2[2]    
            if max(test_redundancy_1[0],test_redundancy_2[0]) >= tolerance: # if the reactions are the same replace add the metabolite to the list of equivalent
                if max(test_redundancy_1[0],test_redundancy_2[0]) >= test_redundancy_ref:
                    test_redundancy_ref = max(test_redundancy_1[0],test_redundancy_2[0])
                    best_match = [xth_n2, y,test_redundancy_1[0], test_redundancy_2[0]]
            if best_match:
                perfect_matches.append(best_match)
                best_match = []


    ''' 3. Change ids and annotation in network 2 to be compatible with network 1 '''

    print('Making annotation compatible between networks')

    equivalent_metabolites_key_list = {x[1]:x[0] for x in n1_n2_equivalent_metabolites} # key:n2 id, value: n1 id
    n2_n1_equivalent_genes = {x[1]:x[0] for x in equivalent_genes} # key:n2 id, value: n1 id

    for reaction_n2 in network_2.reactions:
        # 3.1: Identify metabolites in reaction from network 2 
        metabolites_in_n2 = [x.id for x in reaction_n2.metabolites]
        reaction = reaction_n2.reaction
        for x in metabolites_in_n2:
            x_unique = re.sub(network_2.metabolites.get_by_id(x).compartment+'$', '', network_2.metabolites.get_by_id(x).id)
            if x_unique in equivalent_metabolites_key_list.keys():
                #equivalent_metabolites_key_list.keys()      
                reaction = re.sub(x, equivalent_metabolites_key_list[x_unique]+network_2.metabolites.get_by_id(x).compartment,reaction_n2.reaction)
                reaction_n2.reaction = reaction    
                reaction_n2.metabolites.update()
        # 3.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                reaction_n2.gpr = new_gpr
                reaction_n2.update_genes_from_gpr()

    #network_1_step_3_sec_copy = copy.deepcopy(network_1)
    #network_1 = copy.deepcopy(network_1_step_3_sec_copy)

    ''' 4. Enrich network 1 with the overlaped reactions '''

    # 4.1: Add reaction and metabolite attributes from overlap reaction in netwrok 2 to network 1
    print('Merging overlap reactions')
    
    for x in perfect_matches:
        reaction_n1 = network_1.reactions.get_by_id(x[1])
        reaction_n2 = network_2.reactions.get_by_id(x[0])

        # 2.4.1: Add reaction and metabolite attributes from netwrok 2 to network 1
        reaction = reaction_n2.reaction
        if x[2] < x[3]: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible               
            reaction = reaction.split(' <=> ')[1]+' <=> '+reaction.split(' <=> ')[0]
        if reaction_n1.reversibility == True:
            rev_n2 = re.findall('(...)>', reaction_n2.reaction)[0].strip()+'>'
            reaction = re.sub(rev_n2, '<=>', reaction) #Impose reversibility from network 1

        original_reaction = copy.deepcopy(reaction_n1.reaction)
        reaction_n1.reaction = copy.deepcopy(reaction)
        reaction_n1.metabolites.update()

        # Test the new reaction in the new network

        # 
        #is_imbalanced = get_ids(consistency.find_mass_unbalanced_reactions([reaction_n1]))
        #if is_imbalanced:
        #   eq = reaction_n1.reaction.replace("--", "-").replace("<=", "-") 
        #   mb = RxnBalance2(eq, 'R')
        #    if mb[10] == 0: # re-balanced
        #        species_inverse_lib = {}
        #        for y in species:
        #            if y in Atom_ID_Inverse:
        #                species_inverse_lib[Atom_ID_Inverse[y]] = y # network_1.metabolites.get_by_id(y).id
        #            else:
        #                if y in network_1.metabolites:
        #                    species_inverse_lib[network_1.metabolites.get_by_id(y).formula] = network_1.metabolites.get_by_id(y).id
        #                else:
        #                    species_inverse_lib[network_2.metabolites.get_by_id(y).formula] = network_2.metabolites.get_by_id(y).id

        #        #species_inverse_lib = {network_1.metabolites.get_by_id(x).formula:network_1.metabolites.get_by_id(x).id for x in species}

        #        new_reactnats = mb[6] 
        #        new_products = mb[7]
        #        stch_reactants = mb[0]
        #        stch_products = mb[1]
        #        new_eq = ''
        #        for y in range(0,len(new_reactnats)):
        #            if new_reactnats[y] in species_inverse_lib:
        #                if species_inverse_lib[new_reactnats[y]] in list(network_atom_equivalent.keys()):
        #                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_reactnats[y]]]
        #                else:
        #                    yth_metabolite = species_inverse_lib[new_reactnats[y]]
        #                new_eq = new_eq + str(stch_reactants[y]) + ' ' + yth_metabolite + ' + '
        #            else:
        #                new_eq = new_eq + str(stch_reactants[y]) + ' ' + Atom_ID[new_reactnats[y]] + ' + '
        #        new_eq = re.sub('( [+] $)', ' --> ',new_eq)
        #        for y in range(0,len(new_products)):
        #            if new_products[y] in species_inverse_lib:
        #                if species_inverse_lib[new_products[y]] in list(network_atom_equivalent.keys()):
        #                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_products[y]]]
        #                else:
        #                    yth_metabolite = species_inverse_lib[new_products[y]]
        #                new_eq = new_eq + str(stch_products[y]) + ' ' + yth_metabolite + ' + '
        #            else:
        #                new_eq = new_eq + str(stch_products[y]) + ' ' + Atom_ID[new_products[y]] + ' + '
        #        new_eq = re.sub('( [+]$)', '',new_eq.strip())
        #        x.reaction = new_eq
        #        x.metabolites.update()         
        #model_consistency = consistency.check_stoichiometric_consistency(network_1)
           
        if reaction_n2.id == 'R01009_c': # or is_imbalanced: # model_consistency == False: 
            reaction_n1.reaction = copy.deepcopy(original_reaction)
            reaction_n1.metabolites.update()
            list_of_innconsistent_reactions.append(x[0])
            print('model_inconsistent',x[0])
        else:
            print('model_consistent',x[0])
   
                    
        # 2.4.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))

        reaction_n1._gpr = copy.deepcopy(new_gpr)
        reaction_n1.update_genes_from_gpr()
                            
        # 2.4.2: Add other attributes from netwrok 2 to network 1
        # annotation
        new_annotation = copy.deepcopy(reaction_n1.annotation)
        annotation_overlap = list(reduce(lambda i, j: i & j, (set(x.annotation.keys()) for x in [reaction_n1]+[reaction_n2]))) # common annotations between reaction_n1 and reaction_n2
        for a in annotation_overlap:
            new_annotation[a] = reaction_n2.annotation[a] 
        for a_n2 in reaction_n2.annotation.keys():
            new_annotation[a_n2] = reaction_n2.annotation[a_n2] 

        reaction_n1._annotation = copy.deepcopy(new_annotation)

        #name
        reaction_n1.name = copy.deepcopy(reaction_n2.name)

    #network_1_step_4_sec_copy = copy.deepcopy(network_1)
    #network_2_step_4_sec_copy = copy.deepcopy(network_2)
    #network_1 = copy.deepcopy(network_1_step_4_sec_copy)
    #network_2 = copy.deepcopy(network_2_step_4_sec_copy)
    ''' 5. Enrich network 1 with the non-overlaped reactions '''
    # network_sec_copy_intermediate_5 = copy.deepcopy(network_1)
    reaction_added_from_n2_to_n1 = [x[0] for x in perfect_matches]
    total_reactions_to_test = str(len(network_2.reactions)-len(reaction_added_from_n2_to_n1)-len(imbalanced_reactions))
    for x in network_2.reactions: #[1382:c]:

        # Test mass balance and re-balance if required
        eq = x.reaction  
        eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()
        species = [x.id for x in x.reactants] + [x.id for x in x.products]
        FromList=list()
        for y in species:
            if y in ListOfMetFrom:
                eq = re.sub(y, ListOfMetFrom[y], eq)
            else:
                eq = re.sub(y, Atom_ID_Inverse[y], eq)

        mb = RxnBalance2(eq, 'R') #mb[10] = 0 (re-balanced), 1 (already balanced), 2 (not balanced)

        if mb[10] == 0: # re-balanced
            species_inverse_lib = {}
            for y in species:
                if y in Atom_ID_Inverse:
                    species_inverse_lib[Atom_ID_Inverse[y]] = y # network_1.metabolites.get_by_id(y).id
                else:
                    if y in network_1.metabolites:
                        species_inverse_lib[network_1.metabolites.get_by_id(y).formula] = network_1.metabolites.get_by_id(y).id
                    else:
                        species_inverse_lib[network_2.metabolites.get_by_id(y).formula] = network_2.metabolites.get_by_id(y).id

            #species_inverse_lib = {network_1.metabolites.get_by_id(x).formula:network_1.metabolites.get_by_id(x).id for x in species}

            new_reactnats = mb[6] 
            new_products = mb[7]
            stch_reactants = mb[0]
            stch_products = mb[1]
            new_eq = ''
            for y in range(0,len(new_reactnats)):
                if new_reactnats[y] in species_inverse_lib:
                    if species_inverse_lib[new_reactnats[y]] in list(network_atom_equivalent.keys()):
                        yth_metabolite = network_atom_equivalent[species_inverse_lib[new_reactnats[y]]]
                    else:
                        yth_metabolite = species_inverse_lib[new_reactnats[y]]
                    new_eq = new_eq + str(stch_reactants[y]) + ' ' + yth_metabolite + ' + '
                else:
                    new_eq = new_eq + str(stch_reactants[y]) + ' ' + Atom_ID[new_reactnats[y]] + ' + '
            new_eq = re.sub('( [+] $)', ' --> ',new_eq)
            for y in range(0,len(new_products)):
                if new_products[y] in species_inverse_lib:
                    if species_inverse_lib[new_products[y]] in list(network_atom_equivalent.keys()):
                        yth_metabolite = network_atom_equivalent[species_inverse_lib[new_products[y]]]
                    else:
                        yth_metabolite = species_inverse_lib[new_products[y]]
                    new_eq = new_eq + str(stch_products[y]) + ' ' + yth_metabolite + ' + '
                else:
                    new_eq = new_eq + str(stch_products[y]) + ' ' + Atom_ID[new_products[y]] + ' + '
            new_eq = re.sub('( [+]$)', '',new_eq.strip())
            x.reaction = new_eq
            x.metabolites.update()

        if x not in reaction_added_from_n2_to_n1 and mb[10] != 2: #and not is_imbalanced: # and x.id in imbalanced_reactions:
            try:
                # Give a new id consistent with network1 annotation
                if x.id in network_1.reactions: # in case the id, btt chance, already exist
                    last_reaction = sorted([y.id for y in network_1.reactions],reverse=True)[0]
                    next_reaction = re.findall("^([A-Z,a-z]*?)[0-9]",last_reaction)[0]+str(int((re.findall("^[A-Z,a-z]*?([0-9].*?)$",last_reaction)[0]))+1)
                    x.id = next_reaction

                # Change metabolites id according to network
                metabolites_in_n2 = [y.id for y in x.metabolites]
                reaction = x.reaction
                for y in metabolites_in_n2:
                    # Correct ID in case of miss-annotation
                    if not y in network_1.metabolites and not y in network_2.metabolites:
                        compartment = x.compartments
                        if compartment: compartment = list(compartment)[0]
                        else: compartment = x.id.split('_')[1]
                        new_metabolite_id = y + compartment
                        if not new_metabolite_id in network_1.metabolites: new_metabolite_id = y + '_' + compartment
                        y = new_metabolite_id

                    if not y in network_2.metabolites or not network_2.metabolites.get_by_id(y).compartment:
                        if y in network_1.metabolites:
                            if network_1.metabolites.get_by_id(y).compartment:
                                metabolite_compartment = network_1.metabolites.get_by_id(y).compartment
                            else:
                                metabolite_compartment = list(x.compartments)[0]
                        else:
                            metabolite_compartment = list(x.compartments)[0]
                    else:
                        metabolite_compartment = network_2.metabolites.get_by_id(y).compartment
                    y_unique = re.sub(metabolite_compartment+'$', '', network_2.metabolites.get_by_id(y).id)
                    if y_unique in equivalent_metabolites_key_list.keys():  
                        reaction = re.sub(y, equivalent_metabolites_key_list[y_unique]+network_2.metabolites.get_by_id(y).compartment,x.reaction)
                        x.reaction = reaction    
                        x.metabolites.update()

                # Change gene id according to network_1

                gene_rules = copy.deepcopy(x.gene_reaction_rule)
                gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
                genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
                for g in genes_in_n2:
                    if g in n2_n1_equivalent_genes.keys():                 
                        gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                        new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                        x.gpr = new_gpr
                        x.update_genes_from_gpr()

                # Add the reaction to network_1 and evaluate

                network_1.add_reactions([x])
                #if count < len(network_2.reactions):
                #    model_consistency = True
                #else:
                #    model_consistency = consistency.check_stoichiometric_consistency(network_1)
                model_consistency = consistency.check_stoichiometric_consistency(network_1)
                if model_consistency == False: 
                    network_1.reactions.remove(x.id)
                    list_of_innconsistent_reactions.append(x.id)
                    print('model_inconsistent',x.id)
                else:
                    reaction_added_from_n2_to_n1.append(x.id)
                    print('model_consistent',x.id)
            except:
                print('aborted',x.id)
                continue
        else:
            print('imbalanced',x.id)
        count += 1

    ''' 6. Summary '''

    overlap_reactions = [y[0] for y in perfect_matches]
    not_overlap_reactions = [x.id for x in network_2.reactions if not x.id in [y[0] for y in perfect_matches]]


    return network_1, overlap_reactions, not_overlap_reactions, list_of_innconsistent_reactions


def network_reactions_merge_5 (network_1, network_2, n1_n2_equivalent_metabolites, equivalent_genes): # version 4 + fast addition of metabolic reactions by testing consistency by groups
    """
    Reactions from two metabolic models are compared based on the associated metabolites. If a match is found the information from network_2 replaces the infromation in network_1 or new reactions from network_2 are added to network_1 if no match is found

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    n1_n2_equivalent_metabolites : dictionary with equivalent metabolite ID between network_1 and network_2
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated reaction information from network_2 for the overlaped reactions and new reactions from network_2 in case they were not annotated in network_1 previously
    overlap_reactions: list
    not_overlap_reactions: list
    """

    ''' 1. Warm up'''
    from collections import Counter

    tolerance = 1.75 # values between 0 and 2
    h_in_n_1 = [x.id for x in network_1.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_1
    h_in_n_2 = [x.id for x in network_2.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_2
    ListOfMetFrom = {x.id: x.formula for x in network_1.metabolites}
    ListOfMetFrom_inverse = {x.formula: x.id for x in network_1.metabolites}
      
    Atom_ID = {'H2O1': 'MAM02040', 'H1': 'MAM02039', 'H2O': 'MAM02040', 'H': 'MAM02039', 'Fe': 'MAM01821','F': 'MAM01821', 'X': 'MAM01823', 'R': 'MAM03573', 'Na': 'MAM02519', 'K': 'MAM02200', 'Ca': 'MAM01413'}
    Atom_ID_Inverse = {'MAM02040': 'H2O1', 'MAM02039': 'H1', 'MAM02040': 'H2O', 'MAM02039': 'H', 'MAM01821': 'Fe','MAM01821': 'F', 'MAM01823': 'X', 'MAM03573': 'R', 'MAM02519': 'Na', 'MAM02200': 'K', 'MAM01413': 'Ca'}

    network_atom_equivalent = {'C00000': 'MAM03573'}

    # Evaluate reactions in network_2
    internal_rxns = con_helpers.get_internals(network_2)
    imbalanced_reactions = get_ids(consistency.find_mass_unbalanced_reactions(internal_rxns))
    list_of_innconsistent_reactions = []

    ''' 2. Identify Equivalent Reactions '''

    print('Matching reactions')

    # 2.1. By E.C.
    ec_network_2 = [(x.id, x.annotation['ec-code']) for x in network_2.reactions if 'ec-code' in x.annotation.keys()]
    ec_network_1 = [(x.id, x.annotation['ec-code']) for x in network_1.reactions if 'ec-code' in x.annotation.keys()]
    ec_overlap = list(sorted([[x[0],y[0]] for x in ec_network_2 for y in ec_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.2. By KEGG
    kegg_network_2 = [(x.id, x.annotation['kegg.reaction']) for x in network_2.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_network_1 = [(x.id, x.annotation['kegg.reaction']) for x in network_1.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_overlap = list(sorted([[x[0],y[0]] for x in kegg_network_2 for y in kegg_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.3. Combine KEGG + E.C.
    from more_itertools import locate

    def find_indices(list_to_check, item_to_find):
        indices = locate(list_to_check, lambda x: x == item_to_find)
        return list(indices)

    kegg_ec_overlap = [x for x in kegg_overlap if x in ec_overlap]
    counts = dict(Counter([x[0] for x in kegg_ec_overlap]))
    dm = {key:value for key, value in counts.items() if value > 1}
    uniques_matches = [x for x in kegg_ec_overlap if not x[0] in dm.keys()]
    duplicate_matches = [(x,[kegg_ec_overlap[y][1] for y in find_indices([z[0]for z in kegg_ec_overlap], x)]) for x in dm.keys()]

    # 2.4. Find equivalent reactions between networks

    perfect_matches = []

    for x in uniques_matches+duplicate_matches: # When more than one equivalence from ec+kegg, choose the best

        xth_n2 = x[0]
        xth_n1 = x[1]
        if isinstance(xth_n1,str): xth_n1 = [xth_n1]

        reaction_n2 = network_2.reactions.get_by_id(xth_n2)
        if len(reaction_n2.reactants) > 0:
            substrates_n2 = [x.id for x in reaction_n2.reactants if x.id not in h_in_n_2]
        if len(reaction_n2.products) > 0:
            products_n2 = [x.id for x in reaction_n2.products if x.id not in h_in_n_2]  

        test_redundancy_ref = 0
        best_match = []

        for y in xth_n1:
            reaction_n1 = network_1.reactions.get_by_id(y)
            if len(reaction_n1.reactants) > 0: 
                substrates_n1 = [x.id for x in reaction_n1.reactants if x.id not in h_in_n_1]
            if len(reaction_n1.products) > 0: 
                products_n1 = [x.id for x in reaction_n1.products if x.id not in h_in_n_1]

            if len(reaction_n2.products) == 0 and len(reaction_n1.products) == 0: # reaction in network 1 and 2 are exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]*2
                test_redundancy_2 = list((0,[0],[0]))
            if len(reaction_n2.products) != 0 and len(reaction_n1.products) != 0: # reaction in network 1 and 2 are not exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1) + compare_metabolite_groups(products_n2, products_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]+test_redundancy_1[2]
                del test_redundancy_1[2]
                reversibility_factor = reaction_n2.reversibility*reaction_n1.reversibility # 1: both reactions are reversible, 0: one or both reactions are irreversible
                test_redundancy_2 = list(compare_metabolite_groups(substrates_n2, products_n1, network_2, network_1) + compare_metabolite_groups(products_n2, substrates_n1, network_2, network_1))
                test_redundancy_2[0] = reversibility_factor*(test_redundancy_2[0]+test_redundancy_2[2])
                del test_redundancy_2[2]    
            if max(test_redundancy_1[0],test_redundancy_2[0]) >= tolerance: # if the reactions are the same replace add the metabolite to the list of equivalent
                if max(test_redundancy_1[0],test_redundancy_2[0]) >= test_redundancy_ref:
                    test_redundancy_ref = max(test_redundancy_1[0],test_redundancy_2[0])
                    best_match = [xth_n2, y,test_redundancy_1[0], test_redundancy_2[0]]
            if best_match:
                perfect_matches.append(best_match)
                best_match = []


    ''' 3. Change ids and annotation in network 2 to be compatible with network 1 '''

    print('Making annotation compatible between networks')

    equivalent_metabolites_key_list = {x[1]:x[0] for x in n1_n2_equivalent_metabolites} # key:n2 id, value: n1 id
    n2_n1_equivalent_genes = {x[1]:x[0] for x in equivalent_genes} # key:n2 id, value: n1 id

    for reaction_n2 in network_2.reactions:
        # 3.1: Identify metabolites in reaction from network 2 
        metabolites_in_n2 = [x.id for x in reaction_n2.metabolites]
        reaction = reaction_n2.reaction
        for x in metabolites_in_n2:
            x_unique = re.sub(network_2.metabolites.get_by_id(x).compartment+'$', '', network_2.metabolites.get_by_id(x).id)
            if x_unique in equivalent_metabolites_key_list.keys():
                #equivalent_metabolites_key_list.keys()      
                reaction = re.sub(x, equivalent_metabolites_key_list[x_unique]+network_2.metabolites.get_by_id(x).compartment,reaction_n2.reaction)
                reaction_n2.reaction = reaction    
                reaction_n2.metabolites.update()
        # 3.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                reaction_n2.gpr = new_gpr
                reaction_n2.update_genes_from_gpr()


    ''' 4. Enrich network 1 with the overlaped reactions '''

    # 4.1: Add reaction and metabolite attributes from overlap reaction in netwrok 2 to network 1
    print('Merging overlap reactions')
    
    for x in perfect_matches:
        reaction_n1 = network_1.reactions.get_by_id(x[1])
        reaction_n2 = network_2.reactions.get_by_id(x[0])

        # 2.4.1: Add reaction and metabolite attributes from netwrok 2 to network 1
        reaction = reaction_n2.reaction
        if x[2] < x[3]: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible               
            reaction = reaction.split(' <=> ')[1]+' <=> '+reaction.split(' <=> ')[0]
        if reaction_n1.reversibility == True:
            rev_n2 = re.findall('(...)>', reaction_n2.reaction)[0].strip()+'>'
            reaction = re.sub(rev_n2, '<=>', reaction) #Impose reversibility from network 1

        original_reaction = copy.deepcopy(reaction_n1.reaction)
        reaction_n1.reaction = copy.deepcopy(reaction)
        reaction_n1.metabolites.update()

        # Test the new reaction in the new network

        # 
        #is_imbalanced = get_ids(consistency.find_mass_unbalanced_reactions([reaction_n1]))
        #if is_imbalanced:
        #   eq = reaction_n1.reaction.replace("--", "-").replace("<=", "-") 
        #   mb = RxnBalance2(eq, 'R')
        #    if mb[10] == 0: # re-balanced
        #        species_inverse_lib = {}
        #        for y in species:
        #            if y in Atom_ID_Inverse:
        #                species_inverse_lib[Atom_ID_Inverse[y]] = y # network_1.metabolites.get_by_id(y).id
        #            else:
        #                if y in network_1.metabolites:
        #                    species_inverse_lib[network_1.metabolites.get_by_id(y).formula] = network_1.metabolites.get_by_id(y).id
        #                else:
        #                    species_inverse_lib[network_2.metabolites.get_by_id(y).formula] = network_2.metabolites.get_by_id(y).id

        #        #species_inverse_lib = {network_1.metabolites.get_by_id(x).formula:network_1.metabolites.get_by_id(x).id for x in species}

        #        new_reactnats = mb[6] 
        #        new_products = mb[7]
        #        stch_reactants = mb[0]
        #        stch_products = mb[1]
        #        new_eq = ''
        #        for y in range(0,len(new_reactnats)):
        #            if new_reactnats[y] in species_inverse_lib:
        #                if species_inverse_lib[new_reactnats[y]] in list(network_atom_equivalent.keys()):
        #                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_reactnats[y]]]
        #                else:
        #                    yth_metabolite = species_inverse_lib[new_reactnats[y]]
        #                new_eq = new_eq + str(stch_reactants[y]) + ' ' + yth_metabolite + ' + '
        #            else:
        #                new_eq = new_eq + str(stch_reactants[y]) + ' ' + Atom_ID[new_reactnats[y]] + ' + '
        #        new_eq = re.sub('( [+] $)', ' --> ',new_eq)
        #        for y in range(0,len(new_products)):
        #            if new_products[y] in species_inverse_lib:
        #                if species_inverse_lib[new_products[y]] in list(network_atom_equivalent.keys()):
        #                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_products[y]]]
        #                else:
        #                    yth_metabolite = species_inverse_lib[new_products[y]]
        #                new_eq = new_eq + str(stch_products[y]) + ' ' + yth_metabolite + ' + '
        #            else:
        #                new_eq = new_eq + str(stch_products[y]) + ' ' + Atom_ID[new_products[y]] + ' + '
        #        new_eq = re.sub('( [+]$)', '',new_eq.strip())
        #        x.reaction = new_eq
        #        x.metabolites.update()         
        #model_consistency = consistency.check_stoichiometric_consistency(network_1)
           
        if reaction_n2.id == 'R01009_c': # or is_imbalanced: # model_consistency == False: 
            reaction_n1.reaction = copy.deepcopy(original_reaction)
            reaction_n1.metabolites.update()
            list_of_innconsistent_reactions.append(x[0])
            print('model_inconsistent',x[0])
        else:
            print('model_consistent',x[0])
   
                    
        # 2.4.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))

        reaction_n1._gpr = copy.deepcopy(new_gpr)
        reaction_n1.update_genes_from_gpr()
                            
        # 2.4.2: Add other attributes from netwrok 2 to network 1
        # annotation
        new_annotation = copy.deepcopy(reaction_n1.annotation)
        annotation_overlap = list(reduce(lambda i, j: i & j, (set(x.annotation.keys()) for x in [reaction_n1]+[reaction_n2]))) # common annotations between reaction_n1 and reaction_n2
        for a in annotation_overlap:
            new_annotation[a] = reaction_n2.annotation[a] 
        for a_n2 in reaction_n2.annotation.keys():
            new_annotation[a_n2] = reaction_n2.annotation[a_n2] 

        reaction_n1._annotation = copy.deepcopy(new_annotation)

        #name
        reaction_n1.name = copy.deepcopy(reaction_n2.name)


    ''' 5. Enrich network 1 with the non-overlaped reactions '''

    def reaction_rebalanced (mb, species, Atom_ID_Inverse, network_1, network_atom_equivalent, reaction):
        import numbers
        import decimal
        species_inverse_lib = {}
        for y in species:
            try:
                if isinstance(int(y[-1]), numbers.Number):
                    test_compartment_in_id = False
            except:
                test_compartment_in_id = True
            if test_compartment_in_id == False:
                if y in network_2.metabolites:
                    if network_2.metabolites.get_by_id(y).compartment:
                        y_new = y + network_2.metabolites.get_by_id(y).compartment
                elif y in network_1.metabolites:
                    if network_1.metabolites.get_by_id(y).compartment:
                        y_new = y + network_1.metabolites.get_by_id(y).compartment                   
                else:
                    if y in [x.id for x in x.reactants]:
                        compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.reactants if x.id in network_1.metabolites]))[0]
                        y_new = y + compartment
                    else:
                        compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.products if x.id in network_1.metabolites]))[0]
                        y_new = y + compartment
            else:
                y_new = y

            if y in Atom_ID_Inverse:
                species_inverse_lib[Atom_ID_Inverse[y]] = y_new # network_1.metabolites.get_by_id(y).id
            else:
                if y in network_1.metabolites:
                    species_inverse_lib[network_1.metabolites.get_by_id(y).formula] = network_1.metabolites.get_by_id(y).id
                else:
                    species_inverse_lib[network_2.metabolites.get_by_id(y).formula] = network_2.metabolites.get_by_id(y).id

        #species_inverse_lib = {network_1.metabolites.get_by_id(x).formula:network_1.metabolites.get_by_id(x).id for x in species}

        new_reactnats = mb[6] 
        new_products = mb[7]
        stch_reactants = mb[0]
        stch_products = mb[1]
        new_eq = ''
        for y in range(0,len(new_reactnats)):
            if new_reactnats[y] in species_inverse_lib:
                if species_inverse_lib[new_reactnats[y]] in list(network_atom_equivalent.keys()):
                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_reactnats[y]]]
                else:
                    yth_metabolite = species_inverse_lib[new_reactnats[y]]
                new_eq = new_eq + str(stch_reactants[y]) + ' ' + yth_metabolite + ' + '
            else:
                if reaction.compartments:
                    new_id = Atom_ID[new_reactnats[y]] + str(list(reaction.compartments)[0])
                else:
                    yth_compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.reactants if x.id in network_1.metabolites]))[0]
                    new_id = Atom_ID[new_reactnats[y]] + str(yth_compartment)                    
                new_eq = new_eq + str(stch_reactants[y]) + ' ' + new_id + ' + ' #Atom_ID[new_reactnats[y]] + ' + '

        new_eq = re.sub('( [+] $)', ' --> ',new_eq)
        for y in range(0,len(new_products)):
            if new_products[y] in species_inverse_lib:
                if species_inverse_lib[new_products[y]] in list(network_atom_equivalent.keys()):
                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_products[y]]]
                else:
                    yth_metabolite = species_inverse_lib[new_products[y]]
                new_eq = new_eq + str(stch_products[y]) + ' ' + yth_metabolite + ' + '
            else:
                if reaction.compartments:
                    new_id = Atom_ID[new_products[y]] + str(list(reaction.compartments)[0])
                else:
                    yth_compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.products if x.id in network_1.metabolites]))[0]
                    new_id = Atom_ID[new_products[y]] + str(yth_compartment)
                new_eq = new_eq + str(stch_products[y]) + ' ' + new_id + ' + '

        new_eq = re.sub('( [+]$)', '',new_eq.strip())
        return new_eq


    reaction_added_from_n2_to_n1 = [x[0] for x in perfect_matches]

    count = 0
    batch_size = 200

    pre_batch_network_1 = copy.deepcopy(network_1)
    
    while count < len(network_2.reactions):

        
        if count + batch_size > len(network_2.reactions):  batch_size = len(network_2.reactions) - count # Re-dimension the size of the last batch 

        ''' Testing batchs '''

        batch_pool = [] 
        for x in network_2.reactions[count:count+batch_size]: # Check if the reactions in the batch are mass balanced and add to the batch pool to test in case they are
             
            # Test mass balance and re-balance if required
            eq = x.reaction  
            eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()
            species = [x.id for x in x.reactants] + [x.id for x in x.products]
            for y in species:
                if y in ListOfMetFrom:
                    eq = re.sub(y, ListOfMetFrom[y], eq)
                else:
                    eq = re.sub(y, Atom_ID_Inverse[y], eq)

            mb = RxnBalance2(eq, 'R') #mb[10] = 0 (re-balanced), 1 (already balanced), 2 (not balanced)

            if mb[10] == 0: # re-balanced
                new_eq = reaction_rebalanced (mb, species, Atom_ID_Inverse, network_1, network_atom_equivalent, x)
                x.reaction = new_eq
                x.metabolites.update()

            if mb[10] == 2:
                list_of_innconsistent_reactions.append(x.id)
                print('imbalanced reaction not added to batch pool',x.id)

            if mb[10] != 2: # add reactions to the model and check if it is consistent
                try:
                    # Give a new id consistent with network1 annotation
                    if x.id in network_1.reactions: # in case the id, by chance, already exist in network_1
                        last_reaction = sorted([y.id for y in network_1.reactions],reverse=True)[0]
                        next_reaction = re.findall("^([A-Z,a-z]*?)[0-9]",last_reaction)[0]+str(int((re.findall("^[A-Z,a-z]*?([0-9].*?)$",last_reaction)[0]))+1)
                        x.id = next_reaction

                    # Change metabolites id according to network
                    metabolites_in_n2 = [y.id for y in x.metabolites]
                    reaction = x.reaction
                    for y in metabolites_in_n2:
                        # Correct ID in case of miss-annotation
                        if not y in network_1.metabolites and not y in network_2.metabolites:
                            compartment = x.compartments
                            if compartment: compartment = list(compartment)[0]
                            else: compartment = x.id.split('_')[1]
                            new_metabolite_id = y + compartment
                            if not new_metabolite_id in network_1.metabolites: new_metabolite_id = y + '_' + compartment
                            y = new_metabolite_id

                        if not y in network_2.metabolites or not network_2.metabolites.get_by_id(y).compartment:
                            if y in network_1.metabolites:

                                if network_1.metabolites.get_by_id(y).compartment:
                                    metabolite_compartment = network_1.metabolites.get_by_id(y).compartment
                                else:
                                    metabolite_compartment = list(x.compartments)[0]

                            else:
                                metabolite_compartment = list(x.compartments)[0]

                            network_2.metabolites.get_by_id(y).compartment = metabolite_compartment

                        else:
                            metabolite_compartment = network_2.metabolites.get_by_id(y).compartment

                        y_unique = re.sub(metabolite_compartment+'$', '', network_2.metabolites.get_by_id(y).id)
                        if y_unique in equivalent_metabolites_key_list.keys():
                            # change reaction  
                            reaction = re.sub(y, equivalent_metabolites_key_list[y_unique]+network_2.metabolites.get_by_id(y).compartment,x.reaction)
                            x.reaction = reaction    
                            x.metabolites.update()
                            equivalent_id = equivalent_metabolites_key_list[y_unique]+network_2.metabolites.get_by_id(y).compartment
                        else:
                            equivalent_id = y

                        if not y in network_1.metabolites and not equivalent_id in network_1.metabolites:
                            new_metabolite_to_add_in_network_1 = network_2.metabolites.get_by_id(y).copy()
                            new_metabolite_to_add_in_network_1.id = equivalent_id
                            network_1.add_metabolites(new_metabolite_to_add_in_network_1)

                    # Change gene id according to network_1
                    gene_rules = copy.deepcopy(x.gene_reaction_rule)
                    gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
                    genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
                    for g in genes_in_n2:
                        if g in n2_n1_equivalent_genes.keys():                 
                            gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                            new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                            x.gpr = new_gpr
                            x.update_genes_from_gpr()    
                    
                    # Save the balanced reactions in the batch pool
                    batch_pool = batch_pool + [x]

                    print('reaction added to batch pool',x.id)

                except:
                    print('aborted',x.id)
                    continue

        [pre_batch_network_1.add_reactions([r]) for r in batch_pool] # Add balanced reactions to the pre_batch_network_1

        print(len(pre_batch_network_1.reactions), len(network_1.reactions))

        model_consistency = consistency.check_stoichiometric_consistency(pre_batch_network_1) # Test model consistency

        if model_consistency == True:
            [reaction_added_from_n2_to_n1.append(r.id) for r in batch_pool]
            #network_1 = copy.deepcopy(pre_batch_network_1)
            print('pre_batch_model_consistent')
        else:
            print('pre_batch_model_inconsistent trying individual reactions')
            [pre_batch_network_1.reactions.remove(r.id) for r in batch_pool]
            for r in batch_pool:
                pre_batch_network_1.add_reactions([r])
                model_consistency = consistency.check_stoichiometric_consistency(pre_batch_network_1)

                if model_consistency == False: 
                    pre_batch_network_1.reactions.remove(r.id)
                    list_of_innconsistent_reactions.append(r.id)
                    print('model_inconsistent',r.id)
                else:
                    reaction_added_from_n2_to_n1.append(r.id)
                    print('model_consistent',r.id)

        count = count + batch_size

        print (count, len(network_2.reactions))



    ''' 6. Summary '''

    overlap_reactions = [y[0] for y in perfect_matches]
    not_overlap_reactions = [x.id for x in network_2.reactions if not x.id in [y[0] for y in perfect_matches]]


    return network_1, overlap_reactions, not_overlap_reactions, list_of_innconsistent_reactions


def network_reactions_merge_6 (network_1, network_2, n1_n2_equivalent_metabolites, equivalent_genes): # version 4 + fast addition of metabolic reactions by testing consistency by groups
    """
    Reactions from two metabolic models are compared based on the associated metabolites. If a match is found the information from network_2 replaces the infromation in network_1 or new reactions from network_2 are added to network_1 if no match is found

    Inputs
    ----------
    network_1 and network_2 : metabolic model in Cobra format
    n1_n2_equivalent_metabolites : dictionary with equivalent metabolite ID between network_1 and network_2
    
    Output
    -------
    network_1 : metabolic network in Cobra format. Corresponds to network_1 with updated reaction information from network_2 for the overlaped reactions and new reactions from network_2 in case they were not annotated in network_1 previously
    overlap_reactions: list
    not_overlap_reactions: list
    """

    ''' 1. Warm up'''
    from collections import Counter
    from more_itertools import locate


    def reaction_rebalanced (mb, species, Atom_ID_Inverse, network_1, network_atom_equivalent, reaction):
        import numbers
        import decimal
        species_inverse_lib = {}
        for y in species:
            try:
                if isinstance(int(y[-1]), numbers.Number):
                    test_compartment_in_id = False
            except:
                test_compartment_in_id = True
            if test_compartment_in_id == False:
                if y in network_2.metabolites:
                    if network_2.metabolites.get_by_id(y).compartment:
                        y_new = y + network_2.metabolites.get_by_id(y).compartment
                elif y in network_1.metabolites:
                    if network_1.metabolites.get_by_id(y).compartment:
                        y_new = y + network_1.metabolites.get_by_id(y).compartment                   
                else:
                    if y in [x.id for x in x.reactants]:
                        compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.reactants if x.id in network_1.metabolites]))[0]
                        y_new = y + compartment
                    else:
                        compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.products if x.id in network_1.metabolites]))[0]
                        y_new = y + compartment
            else:
                y_new = y

            if y in Atom_ID_Inverse:
                species_inverse_lib[Atom_ID_Inverse[y]] = y_new # network_1.metabolites.get_by_id(y).id
            else:
                if y in network_1.metabolites:
                    species_inverse_lib[network_1.metabolites.get_by_id(y).formula] = network_1.metabolites.get_by_id(y).id
                else:
                    species_inverse_lib[network_2.metabolites.get_by_id(y).formula] = network_2.metabolites.get_by_id(y).id

        #species_inverse_lib = {network_1.metabolites.get_by_id(x).formula:network_1.metabolites.get_by_id(x).id for x in species}

        new_reactnats = mb[6] 
        new_products = mb[7]
        stch_reactants = mb[0]
        stch_products = mb[1]
        new_eq = ''
        for y in range(0,len(new_reactnats)):
            if new_reactnats[y] in species_inverse_lib:
                if species_inverse_lib[new_reactnats[y]] in list(network_atom_equivalent.keys()):
                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_reactnats[y]]]
                else:
                    yth_metabolite = species_inverse_lib[new_reactnats[y]]
                new_eq = new_eq + str(stch_reactants[y]) + ' ' + yth_metabolite + ' + '
            else:
                if reaction.compartments:
                    new_id = Atom_ID[new_reactnats[y]] + str(list(reaction.compartments)[0])
                else:
                    yth_compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.reactants if x.id in network_1.metabolites]))[0]
                    new_id = Atom_ID[new_reactnats[y]] + str(yth_compartment)                    
                new_eq = new_eq + str(stch_reactants[y]) + ' ' + new_id + ' + ' #Atom_ID[new_reactnats[y]] + ' + '

        new_eq = re.sub('( [+] $)', ' --> ',new_eq)
        for y in range(0,len(new_products)):
            if new_products[y] in species_inverse_lib:
                if species_inverse_lib[new_products[y]] in list(network_atom_equivalent.keys()):
                    yth_metabolite = network_atom_equivalent[species_inverse_lib[new_products[y]]]
                else:
                    yth_metabolite = species_inverse_lib[new_products[y]]
                new_eq = new_eq + str(stch_products[y]) + ' ' + yth_metabolite + ' + '
            else:
                if reaction.compartments:
                    new_id = Atom_ID[new_products[y]] + str(list(reaction.compartments)[0])
                else:
                    yth_compartment = list(set([network_1.metabolites.get_by_id(x.id).compartment for x in reaction.products if x.id in network_1.metabolites]))[0]
                    new_id = Atom_ID[new_products[y]] + str(yth_compartment)
                new_eq = new_eq + str(stch_products[y]) + ' ' + new_id + ' + '

        new_eq = re.sub('( [+]$)', '',new_eq.strip())
        return new_eq


    def find_indices(list_to_check, item_to_find):
        indices = locate(list_to_check, lambda x: x == item_to_find)
        return list(indices)


    tolerance = 1.75 # values between 0 and 2
    h_in_n_1 = [x.id for x in network_1.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_1
    h_in_n_2 = [x.id for x in network_2.metabolites if x.formula == 'H1' or x.formula == 'H' or x.formula == 'H+'] # protons in network_2
    ListOfMetFrom = {x.id: x.formula for x in network_1.metabolites}
    ListOfMetFrom_inverse = {x.formula: x.id for x in network_1.metabolites}
      
    Atom_ID = {'H2O1': 'MAM02040', 'H1': 'MAM02039', 'H2O': 'MAM02040', 'H': 'MAM02039', 'Fe': 'MAM01821','F': 'MAM01821', 'X': 'MAM01823', 'R': 'MAM03573', 'Na': 'MAM02519', 'K': 'MAM02200', 'Ca': 'MAM01413'}
    Atom_ID_Inverse = {'MAM02040': 'H2O1', 'MAM02039': 'H1', 'MAM02040': 'H2O', 'MAM02039': 'H', 'MAM01821': 'Fe','MAM01821': 'F', 'MAM01823': 'X', 'MAM03573': 'R', 'MAM02519': 'Na', 'MAM02200': 'K', 'MAM01413': 'Ca'}

    network_atom_equivalent = {'C00000': 'MAM03573'}

    # Evaluate reactions in network_2
    internal_rxns = con_helpers.get_internals(network_2)
    imbalanced_reactions = get_ids(consistency.find_mass_unbalanced_reactions(internal_rxns))
    list_of_innconsistent_reactions = []

    ''' 2. Identify Equivalent Reactions '''

    print('Matching reactions')

    # 2.1. By E.C.
    ec_network_2 = [(x.id, x.annotation['ec-code']) for x in network_2.reactions if 'ec-code' in x.annotation.keys()]
    ec_network_1 = [(x.id, x.annotation['ec-code']) for x in network_1.reactions if 'ec-code' in x.annotation.keys()]
    ec_overlap = list(sorted([[x[0],y[0]] for x in ec_network_2 for y in ec_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.2. By KEGG
    kegg_network_2 = [(x.id, x.annotation['kegg.reaction']) for x in network_2.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_network_1 = [(x.id, x.annotation['kegg.reaction']) for x in network_1.reactions if 'kegg.reaction' in x.annotation.keys()]
    kegg_overlap = list(sorted([[x[0],y[0]] for x in kegg_network_2 for y in kegg_network_1 if x[1] == y[1] and network_2.reactions.get_by_id(x[0]).compartments == network_1.reactions.get_by_id(y[0]).compartments]))

    # 2.3. Combine KEGG + E.C.

    kegg_ec_overlap = [x for x in kegg_overlap if x in ec_overlap]
    counts = dict(Counter([x[0] for x in kegg_ec_overlap]))
    dm = {key:value for key, value in counts.items() if value > 1}
    uniques_matches = [x for x in kegg_ec_overlap if not x[0] in dm.keys()]
    duplicate_matches = [(x,[kegg_ec_overlap[y][1] for y in find_indices([z[0]for z in kegg_ec_overlap], x)]) for x in dm.keys()]

    # 2.4. Find equivalent reactions between networks

    perfect_matches = []
    for x in uniques_matches+duplicate_matches: # When more than one equivalence from ec+kegg, choose the best

        xth_n2 = x[0]
        xth_n1 = x[1]
        if isinstance(xth_n1,str): xth_n1 = [xth_n1]

        reaction_n2 = network_2.reactions.get_by_id(xth_n2)
        if len(reaction_n2.reactants) > 0:
            substrates_n2 = [x.id for x in reaction_n2.reactants if x.id not in h_in_n_2]
        if len(reaction_n2.products) > 0:
            products_n2 = [x.id for x in reaction_n2.products if x.id not in h_in_n_2]  

        test_redundancy_ref = 0
        best_match = []

        for y in xth_n1:
            reaction_n1 = network_1.reactions.get_by_id(y)
            if len(reaction_n1.reactants) > 0: 
                substrates_n1 = [x.id for x in reaction_n1.reactants if x.id not in h_in_n_1]
            if len(reaction_n1.products) > 0: 
                products_n1 = [x.id for x in reaction_n1.products if x.id not in h_in_n_1]

            if len(reaction_n2.products) == 0 and len(reaction_n1.products) == 0: # reaction in network 1 and 2 are exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]*2
                test_redundancy_2 = list((0,[0],[0]))
            if len(reaction_n2.products) != 0 and len(reaction_n1.products) != 0: # reaction in network 1 and 2 are not exchange reactions
                # Jaccard Index
                test_redundancy_1 = list(compare_metabolite_groups(substrates_n2, substrates_n1, network_2, network_1) + compare_metabolite_groups(products_n2, products_n1, network_2, network_1))
                test_redundancy_1[0] = test_redundancy_1[0]+test_redundancy_1[2]
                del test_redundancy_1[2]
                reversibility_factor = reaction_n2.reversibility*reaction_n1.reversibility # 1: both reactions are reversible, 0: one or both reactions are irreversible
                test_redundancy_2 = list(compare_metabolite_groups(substrates_n2, products_n1, network_2, network_1) + compare_metabolite_groups(products_n2, substrates_n1, network_2, network_1))
                test_redundancy_2[0] = reversibility_factor*(test_redundancy_2[0]+test_redundancy_2[2])
                del test_redundancy_2[2]    
            if max(test_redundancy_1[0],test_redundancy_2[0]) >= tolerance: # if the reactions are the same replace add the metabolite to the list of equivalent
                if max(test_redundancy_1[0],test_redundancy_2[0]) >= test_redundancy_ref:
                    test_redundancy_ref = max(test_redundancy_1[0],test_redundancy_2[0])
                    best_match = [xth_n2, y,test_redundancy_1[0], test_redundancy_2[0]]
            if best_match:
                perfect_matches.append(best_match)
                best_match = []


    ''' 3. Change ids and annotation in network 2 to be compatible with network 1 '''

    print('Making annotation compatible between networks')

    equivalent_metabolites_key_list = {x[1]:x[0] for x in n1_n2_equivalent_metabolites} # key:n2 id, value: n1 id
    n2_n1_equivalent_genes = {x[1]:x[0] for x in equivalent_genes} # key:n2 id, value: n1 id

    for reaction_n2 in network_2.reactions:
        # 3.1: Identify metabolites in reaction from network 2 
        metabolites_in_n2 = [x.id for x in reaction_n2.metabolites]
        reaction = reaction_n2.reaction
        for x in metabolites_in_n2:
            x_unique = re.sub(network_2.metabolites.get_by_id(x).compartment+'$', '', network_2.metabolites.get_by_id(x).id)
            if x_unique in equivalent_metabolites_key_list.keys():
                #equivalent_metabolites_key_list.keys()      
                reaction = re.sub(x, equivalent_metabolites_key_list[x_unique]+network_2.metabolites.get_by_id(x).compartment,reaction_n2.reaction)
                reaction_n2.reaction = reaction    
                reaction_n2.metabolites.update()
        # 3.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                reaction_n2.gpr = new_gpr
                reaction_n2.update_genes_from_gpr()


    ''' 4. Enrich network 1 with the overlaped reactions '''

    # 4.1: Add reaction and metabolite attributes from overlap reaction in netwrok 2 to network 1
    print('Merging overlap reactions')

    #reaction_added_from_n2_to_n1 = [x[0] for x in perfect_matches]
    reaction_added_from_n2_to_n1 = []
    count = 0
    batch_size = 50
  
    sec_copy_network_1 = copy.deepcopy(network_1)

#    network_1 = copy.deepcopy(sec_copy_network_1)

    while count < len(perfect_matches):
        
        if count + batch_size > len(perfect_matches):  batch_size = len(perfect_matches) - count # Re-dimension the size of the last batch 

        ''' Testing batchs '''

        batch_pool = []
        batch_pool_2 = []
 

        for x in perfect_matches[count:count+batch_size]: # Check if the reactions in the batch are mass balanced and add to the batch pool to test in case they are
            if type(x[0]) == list: x = x[0]

            reaction_n2 = network_2.reactions.get_by_id(x[0])
            reaction_n1 = network_1.reactions.get_by_id(x[1])


            # Test the new reaction in the new network
            
            is_imbalanced = get_ids(consistency.find_mass_unbalanced_reactions([reaction_n2]))
            
            if is_imbalanced:
                eq = reaction_n2.reaction  
                eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()
                species = [x.id for x in reaction_n2.reactants] + [x.id for x in reaction_n2.products]
                for y in species:
                    if y in ListOfMetFrom:
                        eq = re.sub(y, ListOfMetFrom[y], eq)
                    else:
                        eq = re.sub(y, Atom_ID_Inverse[y], eq)

                mb = RxnBalance2(eq, 'R') #mb[10] = 0 (re-balanced), 1 (already balanced), 2 (not balanced)

                if mb[10] == 2:
                    list_of_innconsistent_reactions.append(reaction_n2.id)
                    print('imbalanced reaction not added to batch pool',reaction_n2.id)    

                if mb[10] != 2: # add reactions to the model and check if it is consistent    
                    print('reaction added to batch pool',x.id)
                    if mb[10] == 0: # re-balanced
                        new_eq = reaction_rebalanced (mb, species, Atom_ID_Inverse, network_1, network_atom_equivalent, reaction_n2)
                        reaction_n2.reaction = new_eq
                        reaction_n2.metabolites.update()

                    # 2.4.1: Add reaction and metabolite attributes from netwrok 2 to network 1
                    reaction = reaction_n2.reaction
                    if x[2] < x[3]: #check if the reaction is annotated in inverse order in n1 compared with n2 and if both are reversible               
                        reaction = reaction.split(' <=> ')[1]+' <=> '+reaction.split(' <=> ')[0]
                    if reaction_n1.reversibility == True:
                        rev_n2 = re.findall('(...)>', reaction_n2.reaction)[0].strip()+'>'
                        reaction = re.sub(rev_n2, '<=>', reaction) #Impose reversibility from network 1

                    reaction_n2.reaction = copy.deepcopy(reaction)
                    reaction_n2.metabolites.update()

                    #original_reaction = copy.deepcopy(reaction_n1.reaction)
                    reaction_n1.reaction = copy.deepcopy(reaction)
                    reaction_n1.metabolites.update()

                    # Save the balanced reactions in the batch pool
                    batch_pool = batch_pool + [reaction_n1]
                    batch_pool_2 = batch_pool_2 + [reaction_n2]


        model_consistency = consistency.check_stoichiometric_consistency(network_1)

        if model_consistency == True:
            [reaction_added_from_n2_to_n1.append(r.id) for r in batch_pool_2]
            print('pre_batch_model_consistent')

        else:
            print('pre_batch_model_inconsistent trying individual reactions')

            for r in batch_pool:
                network_1.reactions.remove(r.id)
                network_1.add_reactions([sec_copy_network_1.reactions.get_by_id(r.id)])

            for r in batch_pool_2:
                target_reaction_id = [x[1] for x in perfect_matches if x[0] == r.id]
                target_reaction = network_1.reactions.get_by_id(target_reaction_id[0])
                target_reaction.reaction = copy.deepcopy(r.reaction)
                target_reaction.metabolites.update()
                model_consistency = consistency.check_stoichiometric_consistency(network_1)

                if model_consistency == False: 
                    network_1.reactions.remove(r.id)
                    list_of_innconsistent_reactions.append(r.id)
                    print('model_inconsistent',r.id)
                else:
                    reaction_added_from_n2_to_n1.append(r.id)
                    print('model_consistent',r.id)

        count = count + batch_size

    # Update gene annotation
    
    for x in perfect_matches:

        reaction_n1 = network_1.reactions.get_by_id(x[1])
        reaction_n2 = network_2.reactions.get_by_id(x[0])
                  
        # 2.4.2: Add gene attributes from netwrok 2 to network 1              
        # new_gpr = merge_gprs (reaction_n1, reaction_n2) # In case gpr from reaction in n1 and n2 have to be combined uncommend this lines and comment the 3 lines below
        gene_rules = copy.deepcopy(reaction_n2.gene_reaction_rule)
        gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
        genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
        new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
        for g in genes_in_n2:
            if g in n2_n1_equivalent_genes.keys():                 
                gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))

        reaction_n1._gpr = copy.deepcopy(new_gpr)
        reaction_n1.update_genes_from_gpr()
                            
        # 2.4.2: Add other attributes from netwrok 2 to network 1
        # annotation
        new_annotation = copy.deepcopy(reaction_n1.annotation)
        annotation_overlap = list(reduce(lambda i, j: i & j, (set(x.annotation.keys()) for x in [reaction_n1]+[reaction_n2]))) # common annotations between reaction_n1 and reaction_n2
        for a in annotation_overlap:
            new_annotation[a] = reaction_n2.annotation[a] 
        for a_n2 in reaction_n2.annotation.keys():
            new_annotation[a_n2] = reaction_n2.annotation[a_n2] 

        reaction_n1._annotation = copy.deepcopy(new_annotation)

        #name
        reaction_n1.name = copy.deepcopy(reaction_n2.name)


    ''' 5. Enrich network 1 with the non-overlaped reactions '''

    # reaction_added_from_n2_to_n1 = [x[0] for x in perfect_matches]

    count = 0
    batch_size = 200
  
    while count < len(network_2.reactions):
        
        if count + batch_size > len(network_2.reactions):  batch_size = len(network_2.reactions) - count # Re-dimension the size of the last batch 

        ''' Testing batchs '''

        batch_pool = [] 
        for x in network_2.reactions[count:count+batch_size]: # Check if the reactions in the batch are mass balanced and add to the batch pool to test in case they are
             
            # Test mass balance and re-balance if required
            eq = x.reaction  
            eq = re.sub('(^| )[0-9\.]+','', re.sub('<=>', '->', re.sub('-->','->',eq))).strip()
            species = [x.id for x in x.reactants] + [x.id for x in x.products]
            for y in species:
                if y in ListOfMetFrom:
                    eq = re.sub(y, ListOfMetFrom[y], eq)
                else:
                    eq = re.sub(y, Atom_ID_Inverse[y], eq)

            mb = RxnBalance2(eq, 'R') #mb[10] = 0 (re-balanced), 1 (already balanced), 2 (not balanced)

            if mb[10] == 0: # re-balanced
                new_eq = reaction_rebalanced (mb, species, Atom_ID_Inverse, network_1, network_atom_equivalent, x)
                x.reaction = new_eq
                x.metabolites.update()

            if mb[10] == 2:
                list_of_innconsistent_reactions.append(x.id)
                print('imbalanced reaction not added to batch pool',x.id)

            if mb[10] != 2: # add reactions to the model and check if it is consistent
                try:
                    # Give a new id consistent with network1 annotation
                    if x.id in network_1.reactions: # in case the id, by chance, already exist in network_1
                        last_reaction = sorted([y.id for y in network_1.reactions],reverse=True)[0]
                        next_reaction = re.findall("^([A-Z,a-z]*?)[0-9]",last_reaction)[0]+str(int((re.findall("^[A-Z,a-z]*?([0-9].*?)$",last_reaction)[0]))+1)
                        x.id = next_reaction

                    # Change metabolites id according to network
                    metabolites_in_n2 = [y.id for y in x.metabolites]
                    reaction = x.reaction
                    for y in metabolites_in_n2:
                        # Correct ID in case of miss-annotation
                        if not y in network_1.metabolites and not y in network_2.metabolites:
                            compartment = x.compartments
                            if compartment: compartment = list(compartment)[0]
                            else: compartment = x.id.split('_')[1]
                            new_metabolite_id = y + compartment
                            if not new_metabolite_id in network_1.metabolites: new_metabolite_id = y + '_' + compartment
                            y = new_metabolite_id

                        if not y in network_2.metabolites or not network_2.metabolites.get_by_id(y).compartment:
                            if y in network_1.metabolites:

                                if network_1.metabolites.get_by_id(y).compartment:
                                    metabolite_compartment = network_1.metabolites.get_by_id(y).compartment
                                else:
                                    metabolite_compartment = list(x.compartments)[0]

                            else:
                                metabolite_compartment = list(x.compartments)[0]

                            network_2.metabolites.get_by_id(y).compartment = metabolite_compartment

                        else:
                            metabolite_compartment = network_2.metabolites.get_by_id(y).compartment

                        y_unique = re.sub(metabolite_compartment+'$', '', network_2.metabolites.get_by_id(y).id)
                        if y_unique in equivalent_metabolites_key_list.keys():
                            # change reaction  
                            reaction = re.sub(y, equivalent_metabolites_key_list[y_unique]+network_2.metabolites.get_by_id(y).compartment,x.reaction)
                            x.reaction = reaction    
                            x.metabolites.update()
                            equivalent_id = equivalent_metabolites_key_list[y_unique]+network_2.metabolites.get_by_id(y).compartment
                        else:
                            equivalent_id = y

                        if not y in network_1.metabolites and not equivalent_id in network_1.metabolites:
                            new_metabolite_to_add_in_network_1 = network_2.metabolites.get_by_id(y).copy()
                            new_metabolite_to_add_in_network_1.id = equivalent_id
                            network_1.add_metabolites(new_metabolite_to_add_in_network_1)

                    # Change gene id according to network_1
                    gene_rules = copy.deepcopy(x.gene_reaction_rule)
                    gene_rules_2 = gene_rules.replace('and', 'or').replace('(','').replace(')','')
                    genes_in_n2 = list(sorted(gene_rules_2.split(' or ')))
                    for g in genes_in_n2:
                        if g in n2_n1_equivalent_genes.keys():                 
                            gene_rules = re.sub(g, n2_n1_equivalent_genes[g],gene_rules)
                            new_gpr = cobra.core.gene.GPR.from_string(str(gene_rules))
                            x.gpr = new_gpr
                            x.update_genes_from_gpr()    
                    
                    # Save the balanced reactions in the batch pool
                    batch_pool = batch_pool + [x]

                    print('reaction added to batch pool',x.id)

                except:
                    print('aborted',x.id)
                    continue

        [network_1.add_reactions([r]) for r in batch_pool] # Add balanced reactions to the network_1

        print(len(network_1.reactions), len(network_1.reactions))

        model_consistency = consistency.check_stoichiometric_consistency(network_1) # Test model consistency

        if model_consistency == True:
            [reaction_added_from_n2_to_n1.append(r.id) for r in batch_pool]
            #network_1 = copy.deepcopy(network_1)
            print('pre_batch_model_consistent')
        else:
            print('pre_batch_model_inconsistent trying individual reactions')
            [network_1.reactions.remove(r.id) for r in batch_pool]
            for r in batch_pool:
                network_1.add_reactions([r])
                model_consistency = consistency.check_stoichiometric_consistency(network_1)

                if model_consistency == False: 
                    network_1.reactions.remove(r.id)
                    list_of_innconsistent_reactions.append(r.id)
                    print('model_inconsistent',r.id)
                else:
                    reaction_added_from_n2_to_n1.append(r.id)
                    print('model_consistent',r.id)

        count = count + batch_size

        print (count, len(network_2.reactions))



    ''' 6. Summary '''

    overlap_reactions = [y[0] for y in perfect_matches]
    not_overlap_reactions = [x.id for x in network_2.reactions if not x.id in [y[0] for y in perfect_matches]]


    return network_1, overlap_reactions, not_overlap_reactions, list_of_innconsistent_reactions


# Step 4: Delete not used items
'Identify and remove orphan metabolites'
def delete_isolated_metabolites (model):
    """
    Detects metabolites not associated to any reaction (orphan) and removes them

    Inputs
    ----------
    model : metabolic model in Cobra format
    
    Output
    -------
    model : metabolic network in Cobra format with no orphan metabolites
    """
    met_in_rxn = [x.metabolites.keys() for x in model.reactions]
    unique_met_in_rxn = []
    for x in met_in_rxn:
        for y in list(x): unique_met_in_rxn.append(y) 
    unique_met_in_rxn = set(unique_met_in_rxn)
    isolated_metabolites = [x.id for x in model.metabolites if not x in unique_met_in_rxn]
    for x in isolated_metabolites: model.metabolites.get_by_id(x).remove_from_model()
    return model, isolated_metabolites

'Identify and remove reactions with no associated metabolites'
def delete_not_used_reactions (model):
    """
    Detects reactions not associated to any metabolite and removes them

    Inputs
    ----------
    model : metabolic model in Cobra format
    
    Output
    -------
    model : metabolic network in Cobra format without not used reactions
    """
    not_used_reactions = [x.id for x in model.reactions if not x.metabolites]
    model.remove_reactions(not_used_reactions)
    return model

'Identify and remove genes with no associated metabolites'
def delete_not_used_genes (model):
    """
    Detects genes not associated to any reaction and removes them

    Inputs
    ----------
    model : metabolic model in Cobra format
    
    Output
    -------
    model : metabolic network in Cobra format without not used genes
    """  
    unique_gene_in_rxn = [*set([x for x in [x.strip().replace("'","").replace("]","").replace("[","").replace("and","").replace("or","").replace(")","").replace("(","") for x in str([x.gene_reaction_rule.split(" ") for x in model.reactions if x.genes]).split(",")] if x])]
    not_used_genes = [x.id for x in model.genes if x.id not in unique_gene_in_rxn]
    for x in not_used_genes: model.genes.remove(x)
    return model, not_used_genes

