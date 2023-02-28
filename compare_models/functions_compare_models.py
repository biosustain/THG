from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model
import re
from typing import Optional, Union
import xlsxwriter
# for mass balance test
from collections import defaultdict
from functools import reduce


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


def umbalance_test (x, ListOfMetFrom):
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


def defmodel(model,worksheet, workbook, model_x: Optional = None):
    
    row = 1
    bold = workbook.add_format({"bold": True})
    columns = {"A1": "#", "B1": "MODEL","C1":"REFERENCE MODEL"}
    

    for x in columns:
        worksheet.write(x, columns[x], bold)
    

    worksheet.write_column(row, 1, [model.id])
    if model_x != None: worksheet.write_column(row, 2, [model_x.id])


def defgroup(worksheet,workbook):
    row = 1
    bold = workbook.add_format({"bold": True})
    columns = {"A1": "#", "B1": "GROUP", "C1": "ABOUT"}
    about = {1: "no change", 2: "replace", 3: "add", 4:"metabolite in reference model",5:"unique metabolite in reference model",6:"metabolite not in reference model"}
    
    for x in columns:
        worksheet.write(x, columns[x], bold)
        
    for i in range(1,7):
        worksheet.write_column(i, 1, [str(i)])
        worksheet.write_column(i, 1+1, [about[i]])
    

def compartment(model, worksheet, workbook, model_x: Optional = None):
    
    row = 1
    bold = workbook.add_format({"bold": True})
    columns = {"A1": "#", "B1": "COMP", "C1": "COMP ABBR.", "D1": "IN REFERENCE MODEL"}
    
    if model_x == None:
        del columns["D1"]
        
    for x in columns:
        worksheet.write(x, columns[x], bold)
    
    for abbr in model.compartments:
        compartment = model.compartments[abbr]
        worksheet.write_column(row, 1, [compartment])
        worksheet.write_column(row, 2, abbr)
        
        n = "yes"
        if not model_x == None:
            if not abbr in model_x.compartments: 
                n = "no"
            worksheet.write_column(row, 3, [n])
        row += 1


def metabolite(model, worksheet, workbook, model_x: Optional = None):
    row = 0
    MetList = []
    d = {}
    col_num = {"bigg.metabolite": 7, "kegg.compound": 8, "chebi": 9, "vmhmetabolite": 10, "metanetx.chemical": 11, "lipidmaps": 12, "pubchem.compound": 13, "inchikey": 14, "inchi": 15, "hmdb": 16, "sbo":6}
    
    bold = workbook.add_format({"bold": True})
    columns = {"A1": "#", "B1": "ID", "C1": "NAME", "D1": "FORMULA", "E1": "CHARGE", "F1": "COMP ABBR.", "G1": "SBO", "H1": "BIGG", "I1": "KEGG.COMPOUND", "J1": "CHEBI", "K1": "VMH", "L1": "METANETX.CHEMICAL", "M1": "LIPIDMAPS", "N1": "PUBCHEM.COMPOUND", "O1": "INCHIKEY", "P1": "INCHI", "Q1": "HMDB", "R1": "IN REFERENCE MODEL","S1": "REFERENCE MODEL GROUP [4,5,6]"}
    
    if model_x == None:
        del columns["S1"]
        del columns["T1"]
    else: dd = [re.sub(x.compartment, str(), x.id) for x in model_x.metabolites]
        
        
    for x in columns:
        worksheet.write(x, columns[x], bold)
    
    
    for x in model.metabolites:
        x2 = re.sub(x.compartment, str(), x.id)
        #d.setdefault(x2,[]).append(x.compartment)
        
    for x in model.metabolites:
        print(x)
        column = 1
        x2 = re.sub(x.compartment, str(), x.id)
        if x:
            #MetList.append(x2)
            row += 1
            worksheet.write(row, column, x.id)
            #c = ','.join(d[x2])
            
            n = "yes"
            j="4"
            if not model_x == None:
                if not x in model_x.metabolites and x2 in dd:
                    n = x2+" is but not "+x.id
                    #print(x,x2)
                    #print()
                    j="5"
                if not x in model_x.metabolites and not x2 in dd: 
                    n = "no"
                    j="6"
                    #print(x,x2)
                    #print()
                worksheet.write_column(row, 17, [n])
                worksheet.write_column(row, 18, [j])
            
            for y in [x.id, x.name, x.formula, x.charge, x.compartment]:
                worksheet.write(row, column, y)
                column +=1 
            for k, v in x.annotation.items():
                if type(v) == list: v = ",".join(v)
                worksheet.write_column(row, col_num[k], [v])
        

def metabolite_2(model, worksheet, workbook, model_x: Optional = None):
    row = 0
    MetList = []
    d = {}
    col_num = {"bigg.metabolite": 7, "kegg.compound": 8, "chebi.compound": 9, "vmhmetabolite": 10, "metanetx.chemical": 11, "lipidmaps": 12, "pubchem.compound": 13, "inchikey": 14, "inchi": 15, "hmdb": 16, "lipidbank": 17, "sbo":6}
    
    bold = workbook.add_format({"bold": True})
    columns = {"A1": "#", "B1": "ID", "C1": "NAME", "D1": "FORMULA", "E1": "CHARGE", "F1": "COMP ABBR.", "G1": "SBO", "H1": "BIGG", "I1": "KEGG.COMPOUND", "J1": "CHEBI.COMPOUND", "K1": "VMH", "L1": "METANETX.CHEMICAL", "M1": "LIPIDMAPS", "N1": "PUBCHEM.COMPOUND", "O1": "INCHIKEY", "P1": "INCHI", "Q1": "HMDB", "R1": "LIPIDBANK", "S1": "IN REFERENCE MODEL","T1": "REFERENCE MODEL GROUP [4,5,6]"}
    
    if model_x == None:
        del columns["T1"]
        del columns["U1"]
    else: dd = [re.sub(x.compartment, str(), x.id) for x in model_x.metabolites]
        
        
    for x in columns:
        worksheet.write(x, columns[x], bold)
    
    
    for x in model.metabolites:
        x2 = re.sub(x.compartment, str(), x.id)
        #d.setdefault(x2,[]).append(x.compartment)
        
    for x in model.metabolites:
        column = 1
        x2 = re.sub(x.compartment, str(), x.id)
        if x:
            #MetList.append(x2)
            row += 1
            worksheet.write(row, column, x.id)
            #c = ','.join(d[x2])
            
            n = "yes"
            j="4"
            if not model_x == None:
                if not x in model_x.metabolites and x2 in dd:
                    n = x2+" is but not "+x.id
                    #print(x,x2)
                    #print()
                    j="5"
                if not x in model_x.metabolites and not x2 in dd: 
                    n = "no"
                    j="6"
                    #print(x,x2)
                    #print()
                worksheet.write_column(row, 18, [n])
                worksheet.write_column(row, 19, [j])
            
            for y in [x.id, x.name, x.formula, x.charge, x.compartment]:
                worksheet.write(row, column, y)
                column +=1 
            for k, v in x.annotation.items():
                if type(v) == list: v = ",".join(v)
                worksheet.write_column(row, col_num[k], [v])
        

def reaction(model, worksheet, workbook, model_x: Optional = None):
    ListOfMetFrom = {x.id: x.formula for x in model.metabolites}
    row = 0
    col = 0
    bold = workbook.add_format({"bold": True})
    columns = {"A1": "#", "B1": "ID", "C1": "EQUATION_FORMULA", "D1":'EQUATION', "E1": "EC-NUMBER", "F1": "KEGG.REACTION", "G1": "LOWER_BOUND", "H1": "UPPER_BOUND", "I1": 'GENE ASSOCIATION', "J1": "COMP ABBR.", "K1": "IN REFERENCE MODEL", "L1": "MASS_BALANCE GROUP [1,2]", "M1": "KEGG.REACTION GROUP [1,2,3]", "N1": "GENE ASSOCIATION GROUP [1,2,3]"}
    
    if model_x == None:
        del columns["K1"]
        del columns["L1"]
        del columns["M1"]
        del columns["N1"]
        
    for x in columns:
        worksheet.write(x, columns[x], bold)
    
    for x in model.reactions:
        x2 = x.id
        row += 1
        reaction2 = x.reaction
        eq = reaction2
        for y in x.metabolites:
            if ListOfMetFrom[y.id]:
                eq = re.sub(y.id, ListOfMetFrom.get(y.id, y.id), eq)
        
        compartment = ",".join(list(set([y.compartment for y in x.metabolites])))
        z = x.annotation.get('kegg.reaction', '')   
        gpr = x.gene_reaction_rule
        
        n = "yes"
        n2 = "1" # no change stch 
        n3 = "0"
        n4 = "0"
        if not model_x == None:
            if not x in model_x.reactions: 
                n = "no" # new rxn yes "not in ref model"
            else:
                if [*x.metabolites.values()] != [*model_x.reactions.get_by_id(x.id).metabolites.values()] and not umbalance_test (model_x.reactions.get_by_id(x.id), ListOfMetFrom)[0]:
                    n2 = "2" # stch replaced
                if z and z == model_x.reactions.get_by_id(x.id).annotation.get('kegg.reaction', ''):
                    n3 = "1" # kegg no change
                if z and z != model_x.reactions.get_by_id(x.id).annotation.get('kegg.reaction', ''):
                    n3 = "2" # kegg replace
                if not z and model_x.reactions.get_by_id(x.id).annotation.get('kegg.reaction', ''):
                    n3 = "3" # kegg add
                if gpr and gpr == model_x.reactions.get_by_id(x.id).gene_reaction_rule:
                    n4 = "1" # gpr no change
                if gpr and gpr != model_x.reactions.get_by_id(x.id).gene_reaction_rule:
                    n4 = "2" # gpr replaxe 
                if not gpr and model_x.reactions.get_by_id(x.id).gene_reaction_rule:
                    n4 = "3" # gpr add
            worksheet.write_column(row, 10, [n])
            worksheet.write_column(row, 11, n2)
            worksheet.write_column(row, 12, n3)
            worksheet.write_column(row, 13, n4)
            
        y = x.annotation.get('ec-code', '')
        col = 1
        for k in [x2, eq, reaction2, y, z, x.lower_bound, x.upper_bound,gpr, compartment]:
            if type(k) == list: k = ','.join(k)
            worksheet.write_column(row, col, [k])
            col += 1


def gene(model, worksheet, workbook, model_x: Optional = None):
    
    
    row = 1
    column = 1   
    bold = workbook.add_format({"bold": True})
    columns = {"A1": "#", "B1": "NAME", "C1": "ENSEMBLE", "D1":'UNIPROT', "E1": "HGNC.SYMBOL", "F1": "NCBIGENE", "G1": "SBO", "H1": "IN REFERENCE MODEL", "I1": "ENSEMBLE GROUP [1,2,3]"}
    col = {"sbo": 6, "ensembl": 2, "uniprot": 3, "hgnc.symbol": 4, "ncbigene": 5} 
    if model_x == None:
        del columns["H1"]
        del columns["I1"]
        
    for x in columns:
        worksheet.write(x, columns[x], bold)
    
    for x in model.genes:
        g = x.annotation.get('ensembl', '') 
        worksheet.write(row, column, x.id)
        
        n = "yes"
        n2 = "0"
        if not model_x == None:
            if not x in model_x.genes: 
                n = "no"
            worksheet.write_column(row, 7, [n])
            if g and g == model_x.genes.get_by_id(x.id).annotation.get('ensembl', ''): # No change 
                n2 = "1"
            if g and g != model_x.genes.get_by_id(x.id).annotation.get('ensembl', ''): # replace
                n2 = "2"
            if not g and model_x.genes.get_by_id(x.id).annotation.get('ensembl', ''): # add
                n2 = "3"
            worksheet.write_column(row, 8, n2)
                
            
        for k, v in x.annotation.items():
            if type(v) == list: v = ",".join(v)
            worksheet.write_column(row, col[k], [v])
        row += 1
        

def metabolite2(model2, model, worksheet,workbook): # the first input is actually "model" 
    row = 0 
    column = 1
    MetList = []
    col_num = {"bigg.metabolite": 4, "kegg.compound": 3, "chebi": 5, "vmhmetabolite": 6, "metanetx.chemical": 7, "lipidmaps":8, "pubchem.compound": 9, "inchikey": 10, "inchi": 11, "hmdb": 12, "sbo":3}
        
    col2 = {"bigg.metabolite": 'BIGG', "kegg.compound": 'KEGG.COMPOUND', "chebi": 'CHEBI', "vmhmetabolite": 'VMH', "metanetx.chemical": 'METANETX.CHEMICAL', "lipidmaps":'LIPIDMAPS', "pubchem.compound": 'PUBCHEM.COMPOUND', "inchikey": 'INCHIKEY', "inchi": 'INCHIKEY', "hmdb": 'HMDB', "sbo":'SBO'}
    #columns = {"A1": "#", "B1": "ID", "C1": "NAME", "D1": "FORMULA", "E1": "CHARGE", "F1": "COMP ABBR", "G1": "SBO", "H1": "BIGG", "I1": "KEGG.COMPOUND", "J1": "CHEBI", "K1": "VMH", "L1": "METANETX.CHEMICAL", "M1": "LIPIDMAPS", "N1": "PUBCHEM.COMPOUND", "O1": "INCHIKEY", "P1": "INCHI", "Q1": "HMDB", "R1": "SBO", "S1": "IN REFERENCE MODEL","T1": "REFERENCE MODEL GROUP [4,5,6]"}
    columns = {"A1": "#", "B1": "ID", "C1": "FORMULA"}
    bold = workbook.add_format({"bold": True})  
    
    for x in columns:
        worksheet.write(x, columns[x], bold)
    for x in model2.metabolites:
        x2 = re.sub(x.compartment, str(), x.id)
        if not x2 in MetList:
            MetList.append(x2)
            row += 1
            worksheet.write(row, column, x2)
            if x in model.metabolites:
                j="1"
                if x.formula != model.metabolites.get_by_id(x.id).formula:
                    j="2"
                worksheet.write(row, 2, j)
                for y in x.annotation.items():
                    #print(y)
                    num = "0" # no id
                    worksheet.write_column(row, col_num[y[0]], num)
                    if y[0] in col_num:
                        worksheet.write(0, col_num[y[0]], col2[y[0]],bold)
                        #print(y[0],model.metabolites.get_by_id(x.id).annotation)
                        if not y[0] in model.metabolites.get_by_id(x.id).annotation: num = '3'  # id added 
                        elif y[0] in model.metabolites.get_by_id(x.id).annotation and not y[1] == model.metabolites.get_by_id(x.id).annotation[y[0]]: 
                            num = '2' # id replaced 
                        elif y[0] in model.metabolites.get_by_id(x.id).annotation and y[1] == model.metabolites.get_by_id(x.id).annotation[y[0]]: num = '1' # id no change 
                        worksheet.write_column(row, col_num[y[0]], num)


def metabolite2_2(model2, model, worksheet,workbook): # the first input is actually "model" 
    row = 0 
    column = 1
    MetList = []
    col_num = {"bigg.metabolite": 3, "kegg.compound": 4, "chebi.compound": 5, "vmhmetabolite": 6, "metanetx.chemical": 7, "lipidmaps":8, "pubchem.compound": 9, "inchikey": 10, "inchi": 11, "hmdb": 12, "lipidbank": 13, "sbo":14}
        
    col2 = {"bigg.metabolite": 'BIGG', "kegg.compound": 'KEGG.COMPOUND', "chebi.compound": 'CHEBI.COMPOUND', "vmhmetabolite": 'VMH', "metanetx.chemical": 'METANETX.CHEMICAL', "lipidmaps":'LIPIDMAPS', "pubchem.compound": 'PUBCHEM.COMPOUND', "inchikey": 'INCHIKEY', "inchi": 'INCHIKEY', "hmdb": 'HMDB', "lipidbank": 'LIPIDBANK', "sbo":'SBO'}
    columns = {"A1": "#", "B1": "ID", "C1": "FORMULA"}
    bold = workbook.add_format({"bold": True})  
    for x in columns:
        worksheet.write(x, columns[x], bold)
    for x in model2.metabolites:
        x2 = re.sub(x.compartment, str(), x.id)
        if not x2 in MetList:
            MetList.append(x2)
            row += 1
            worksheet.write(row, column, x2)
            if x in model.metabolites:
                j="1"
                if x.formula != model.metabolites.get_by_id(x.id).formula:
                    j="2"
                worksheet.write(row, 2, j)
                for y in x.annotation.items():
                    num = "0" # no id
                    worksheet.write_column(row, col_num[y[0]], num)
                    if y[0] in col_num:
                        worksheet.write(0, col_num[y[0]], col2[y[0]],bold)
                        if not y[0] in model.metabolites.get_by_id(x.id).annotation: num = '3'  # id added 
                        elif y[0] in model.metabolites.get_by_id(x.id).annotation and not y[1] == model.metabolites.get_by_id(x.id).annotation[y[0]]: 
                            num = '2' # id replaced 
                        elif y[0] in model.metabolites.get_by_id(x.id).annotation and y[1] == model.metabolites.get_by_id(x.id).annotation[y[0]]: num = '1' # id no change 
                        worksheet.write_column(row, col_num[y[0]], num)
            else:
                j = "2"
                worksheet.write(row, 2, j)
                for y in x.annotation.items():
                    num = "0" # no id
                    worksheet.write_column(row, col_num[y[0]], num)
                    if y[0] in col_num:
                        worksheet.write(0, col_num[y[0]], col2[y[0]],bold)
                        num = '3'  # id added 
                        worksheet.write_column(row, col_num[y[0]], num)


