#!/usr/bin/python
# -*- coding: utf-8 -*-
#import urllib2
import re
#import numpy as np
#from scipy import linalg
import time
import copy
import pubchempy as pcp
from class_generate_database import *
from pattern_generate_database import *
from equations_generate_database import *
import urllib.request, urllib.error, urllib.parse#,cookielib
import re
import urllib.request, urllib.parse, urllib.error
import requests
import pickle
#from Class import gpr


def mar_number(kegg_id):
    global mar_number_dict
    global mar_number_i
    if not 'mar_number_dict' in globals() and not 'mar_number_i' in globals():
        mar_number_dict = {}
        mar_number_i = 0
    if kegg_id in mar_number_dict:
        mar_number_i = mar_number_dict[kegg_id]
    else:
        mar_number_i = mar_number_i + 1
        mar_number_dict[kegg_id] = mar_number_i
    return 'MAR'+str(mar_number_dict[kegg_id])

def create_dict():
    Compartment_CLpkl = "pkl/Compvariable"
    LocVar = open(Compartment_CLpkl, "rb")
    LocVar = pickle.load(LocVar)
    return LocVar

def create_dict2(Gene,ENSG):
    n=''
    global my_dict
    if 'my_dict' not in globals():
        my_dict = {}
    if Gene and ENSG:
        my_dict[Gene]=ENSG 
    if Gene and not ENSG:
        n= my_dict.get(Gene)
    return my_dict,n


def mam_number(kegg_id):
    global mam_number_dict
    global mam_number_i
    if not 'mam_number_dict' in globals() and not 'mam_number_i' in globals():
        mam_number_dict = {}
        mam_number_i = 0
    if kegg_id in mam_number_dict:
        mam_number_i = mam_number_dict[kegg_id]
    else:
        mam_number_i = mam_number_i + 1
        mam_number_dict[kegg_id] = mam_number_i
    return str(mam_number_dict[kegg_id])


############### Head ############### 
def Head(ModID,ModName):
	head = '''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" xmlns:groups="http://www.sbml.org/sbml/level3/version1/groups/version1" level="3" version="1" fbc:required="false" groups:required="false">
  <model metaid="HumanGEM" id="HumanGEM" name="blankName" fbc:strict="false">
    <notes>
      <body xmlns="http://www.w3.org/1999/xhtml">
        <p>Genome-scale metabolic models are valuable tools to study metabolism and provide a scaffold for the integrative analysis of omics data. This is the latest version of Human-GEM, which is a genome-scale metabolic model of a generic human cell. The objective of Human-GEM is to serve as a community model for enabling integrative and mechanistic studies of human metabolism.</p>
      </body>
    </notes>
    <annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:vCard4="http://www.w3.org/2006/vcard/ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
        <rdf:Description rdf:about="#HumanGEM">
          <dcterms:creator>
            <rdf:Bag/>
          </dcterms:creator>
          <dcterms:created rdf:parseType="Resource">
            <dcterms:W3CDTF>'''+time.strftime("%Y-%m-%dT%H:%M:%SZ")+'''</dcterms:W3CDTF>
          </dcterms:created>
          <dcterms:modified rdf:parseType="Resource">
            <dcterms:W3CDTF>'''+time.strftime("%Y-%m-%dT%H:%M:%SZ")+'''</dcterms:W3CDTF>
          </dcterms:modified>
          <bqbiol:is>
            <rdf:Bag>
              <rdf:li rdf:resource="https://identifiers.org/taxonomy/9606"/>
            </rdf:Bag>
          </bqbiol:is>
        </rdf:Description>
      </rdf:RDF>
    </annotation>
    <listOfUnitDefinitions>
      <unitDefinition id="mmol_per_gDW_per_hr">
        <listOfUnits>
          <unit kind="mole" exponent="1" scale="-3" multiplier="1"/>
          <unit kind="gram" exponent="-1" scale="0" multiplier="1"/>
          <unit kind="second" exponent="-1" scale="0" multiplier="3600"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>'''
	return head



############### Subcel. Loc ############### 
def DefCompart(CSL): 
#	if len(re.sub(' $','',re.sub('^ ','',CSL)).split(" "))>1: 
#		ID = (CSL.split(" ")[0][0]+CSL.split(" ")[1][0]).lower().replace(" ","")       
#	else:
#		if len(CSL.split(" "))>1: 
#			ID = CSL[0:3].lower().replace(" ","")
#		else:
#			ID = CSL[0:2].lower().replace(" ","")
	#LocDict = {"Extracellular" : "s","Peroxisome" : "p","Mitochondria" : "m","Cytosol" : "c","Lysosome" : "l","Endoplasmic reticulum" : "r", "Golgi apparatus" : "g","Nucleus":"n","Inner mitochondria":"i"} 
	#ID = LocDict.get(CSL)   
#	if len(CSL.split(" "))>1: 
#		ID = (CSL.split(" ")[0][0]+CSL.split(" ")[1][0]).lower().replace(" ","")
#	else:
#		ID = CSL[0:2].lower().replace(" ","")
#	ID = CSL[0:2].lower()
	LocVar = create_dict()
	ID = LocVar.get(CSL)
	GO = '0005576'
	CL = '''      <compartment metaid="'''+ID+'''" sboTerm="SBO:0000290" id="'''+ID+'''" name="'''+CSL[0].upper()+CSL[1:]+'''" spatialDimensions="3" size="1" constant="true"/>'''
# GO not in human1 
	return CL


############### Compound ############### 
def DefComp(Comp):    
	CompID = Comp.ID1	  
	#print(136, CompID)    
	#CompIDRef = '<rdf:li rdf:resource="urn:miriam:kegg.compound:'+CompID.split('_')[0]+'" />'
	CompName = Comp.Name()
	if CompName: CompName = CompName.replace("&","").replace(";","")
	#if not CompName:
	#	print(CompID)
	#	#print(0, ID)                                                  
	#LocDict = {"Extracellular" : "s","Peroxisome" : "p","Mitochondria" : "m","Cytosol" : "c","Lysosome" : "l","Endoplasmic reticulum" : "r", "Golgi apparatus" : "g","Nucleus":"n","Inner mitochondria":"i"}
	#print(Comp.Subcel)    
#	Comparment = LocDict.get(Comp.Subcel)
#	Comparment = Comp.Subcel[0:2].lower()
	#if len(re.sub(' $','',re.sub('^ ','',Comp.Subcel)).split(" "))>1: 
	#	Comparment = (Comp.Subcel.split(" ")[0][0]+Comp.Subcel.split(" ")[1][0]).lower().replace(" ","")       
	#else:
	#	if len(Comp.Subcel.split(" "))>1: 
	#		Comparment = Comp.Subcel[0:3].lower().replace(" ","")
	#	else:
	#		Comparment = Comp.Subcel[0:2].lower().replace(" ","")
	#if Comp.Subcel == "Cytosol": Comparment = 'c'
	#if Comp.Subcel == "Extracellular": Comparment = 's'
	#if Comp.Subcel == "Peroxisome": Comparment = 'p'
	#if Comp.Subcel == "Mitochondria": Comparment = 'm'
	#if Comp.Subcel == "Lysosome": Comparment = 'l'
	#if Comp.Subcel == "Endoplasmic reticulum": Comparment = 'r'
	#if Comp.Subcel == "Golgi apparatus": Comparment = 'g'
	#if Comp.Subcel == "Nucleus": Comparment = 'n'
	#if Comp.Subcel == "Inner mitochondria"Comparment = 'i'
	LocVar = create_dict()
	Comparment = LocVar.get(Comp.Subcel[0].upper()+Comp.Subcel[1:])
	Formula = Comp.Formula1().replace('-','') # in case formula from pubchem   
	if "(" in Formula: # in case glycan 
	#	Path="Formula"        
	#	with open(Path,'a') as Formula_file:
	#		Formula_file.write(CompID.split("_")[0]+"\t"+Formula+"\n")
		Formula=Comp.Formula4()       
	charge = Comp.charge()           
	PubChem = Comp.PubChem() # NO
	if PubChem:
		PubChemID = '                  <rdf:li rdf:resource="https://identifiers.org/pubchem.substance/'+PubChem+'"/>\n'         
	else:
		PubChemID = ''         
	CheBI = Comp.CheBI()
	if CheBI:
		CheBIID = '                  <rdf:li rdf:resource="https://identifiers.org/chebi/CHEBI:'+CheBI+'"/>\n'
	else:
		CheBIID = ''
	LIPIDMAPS =  Comp.LIPIDMAPS()
	if LIPIDMAPS:
		LIPIDMAPSID = '                  <rdf:li rdf:resource="https://identifiers.org/lipidmaps/'+LIPIDMAPS+'"/>\n'
	else:
		LIPIDMAPSID = ''
	LipidBank =  Comp.LipidBank() # NO
	if LipidBank:
		LipidBankID = '                  <rdf:li rdf:resource="https://identifiers.org/LipidBank:LipidBank:'+LipidBank+'"/>'
	else:
		LipidBankID = ''
	if Comp.ID1[0] == 'C':
		Kegg = '<rdf:li rdf:resource="https://identifiers.org/kegg.compound/'+re.findall('C[0-9]+',Comp.ID1)[0]+'"/>\n'
	elif Comp.ID1[0] == 'G':
		Kegg = '<rdf:li rdf:resource="https://identifiers.org/kegg.compound/'+re.findall('G[0-9]+',Comp.ID1)[0]+'"/>\n'
	else:
		Kegg = ''       
	if Comp.CID():        
		PubChem_CID = '                  <rdf:li rdf:resource="https://identifiers.org/pubchem.compound/'+Comp.CID()+'"/>\n'
	else:
		PubChem_CID =''        
	if Comp.inchi():
		inchi = '                  <rdf:li rdf:resource="https://identifiers.org/inchi/'+Comp.inchi()+'"/>\n'
	else:
		inchi=''        
	if Comp.inchikey():        
		inchikey = '                  <rdf:li rdf:resource="https://identifiers.org/inchikey/'+Comp.inchikey()+'"/>\n'  
	else:
		inchikey=''                   
	ID = Kegg + PubChemID + PubChem_CID + inchi + inchikey + CheBIID + LIPIDMAPSID + LipidBankID
	KeggComp = Kegg+Comparment    
	ID = ID#[:-2]
	no = mam_number(CompID.split('_')[0]) 
	ModelComp ='''      <species metaid="MAM'''+no+Comparment+'''" sboTerm="SBO:0000247" id="MAM'''+no+Comparment+'''" name="'''+CompName+'''" compartment="'''+Comparment+'''" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false" fbc:charge="'''+charge+'''" fbc:chemicalFormula="'''+Formula+'''">
        <annotation>
          <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:vCard4="http://www.w3.org/2006/vcard/ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
            <rdf:Description rdf:about="#MAM'''+no+Comparment+'''">
              <bqbiol:is>
                <rdf:Bag>
                  '''+ID.strip()+'''
                </rdf:Bag>
              </bqbiol:is>
            </rdf:Description>
          </rdf:RDF>
        </annotation>
      </species>'''
	return ModelComp

def DefDefRxn(Rxn,MetEquiv):
	RxnID = Rxn.ID
	return RxnID


def listOfParameters():
	parameter = '''    <listOfParameters>
      <parameter id="FB1N1000" value="-1000" constant="true"/>
      <parameter id="FB2N0" value="0" constant="true"/>
      <parameter id="FB3N1000" value="1000" constant="true"/>
    </listOfParameters>'''
	return parameter     

############### Reaction ############### 
def DefRxn(Rxn,MetEquiv): 
	try:    
		LocVar = create_dict()
		#print()    
		#print(10000000, Rxn, Rxn.GPR, 0)
		#print()     
		RxnID = Rxn.ID 
		print(RxnID)
		CL =  re.findall('_([A-Za-z ]+)', RxnID)[0]
		#print(CL)    
		RxnName = Rxn.Name()#.replace(";","")
		if RxnName:
			RxnName = Rxn.Name()[0]#.replace(";","")
		if not RxnName: 
			RxnName = RxnID
		else:
			if len(RxnName) < 2 and Rxn.Name()[1]:
				RxnName = Rxn.Name()[1]#.replace(";","")
			else:
				RxnName = Rxn.Name()[0]#.replace(";","")
#			RxnName = Rxn.Name()[0]#.replace(";","")
		#Comparment = Rxn.Subcel[0:2].lower()
		RxnEC = Rxn.EC()
		#print('ø', RxnEC,RxnEC[0])        
		#print(RxnEC)    
		if RxnEC[0]:
			RxnEC = RxnEC[0]
			#print(RxnEC)            
			RxnEC = '''<rdf:li rdf:resource="https://identifiers.org/ec-code/'''+RxnEC+'''"/>'''        
			#RxnEC = RxnEC.replace('+','')
			#print(RxnEC)        
		else: 
			RxnEC = ''
		RxnTermodyn = Rxn.Termodyn()
		#print(RxnID)     
		if RxnTermodyn == 'irreversible':
			RxnRev = 'false'
			lb = 'FB2N0'
		else:
			RxnRev = 'true'
			lb = 'FB1N1000'
		#print(RxnTermodyn, lb, RxnID)    
		RxnPathway = Rxn.Pathway()
		S = []
		corrected_susbtrate = []  #coorective patch to have the rght substrates IDs in RxnList_CL
		for x in Rxn.Substrate:
#			            
#			if not "_" in x[2]: print(RxnID)            
#			if not "_" in x[2]: print(x[2])
			x[2] = x[1].split(':')[-1]+'_'+str(RxnID.split('_')[1])[0].upper()+str(RxnID.split('_')[1])[1:]
			print(x)            
			met = x[2].split("_")[0]        
			cl = x[2].split("_")[1]
			if met in MetEquiv:
				corrected_x2 = str(MetEquiv[met])+"_"+str(cl)
				S = S + [[corrected_x2,str(x[0])]]
			else:
				corrected_x2 = str(met)+"_"+str(cl)
				S = S + [[corrected_x2,str(x[0])]]
			corrected_x = copy.deepcopy(x)
			corrected_x[2] = corrected_x2
			corrected_susbtrate.append(corrected_x) 
		Rxn.SetSubstrate(corrected_susbtrate)
#		RxnSubstrate =  ['<speciesReference species="'+x[0]+'" stoichiometry="'+str(x[1])+'" constant="true"/>' for x in S]
		RxnSubstrate =  ['          <speciesReference species="'+re.sub('_.+','_'+LocVar.get(re.search('_(.+)',x[0]).group(1)),x[0])+'" stoichiometry="'+str(x[1])+'" constant="true"/>\n' for x in S]# ID + Stoichometry
#	RxnSubstrate =  ['<speciesReference species="M_'+x[2]+'" stoichiometry="'+str(x[0])+'"/>' for x in Rxn.Substrate] # ID + Stoichometry
		SubsRef = ''
		z = 0
		while z < len(RxnSubstrate):
			SubsRef = SubsRef + RxnSubstrate[z]
			z = z+1       
		P = []
		corrected_product = []  #coorective patch to have the rght substrates IDs in RxnList_CL
		for x in Rxn.Product:
#			if not "_" in x[2]: print(RxnName)            
#			if not "_" in x[2]: print(RxnID)            
#			if not "_" in x[2]: print(x[2])
			x[2] = x[1].split(':')[-1]+'_'+str(RxnID.split('_')[1])[0].upper()+str(RxnID.split('_')[1])[1:]   
			print(x)    
			met = x[2].split("_")[0]
			cl = x[2].split("_")[1]
			if met in MetEquiv:
				corrected_x2 = str(MetEquiv[met])+"_"+str(cl)
				P = P + [[corrected_x2,str(x[0])]]
			else:
				corrected_x2 = str(met)+"_"+str(cl)
				P = P + [[corrected_x2,str(x[0])]]
			corrected_x = copy.deepcopy(x)
			corrected_x[2] = corrected_x2
			corrected_product.append(corrected_x)
		Rxn.SetProduct(corrected_product)
#		RxnProduct =  ['<speciesReference species="'+x[0]+'" stoichiometry="'+str(x[1])+'" constant="true"/>' for x in P]
		RxnProduct =  ['          <speciesReference species="'+re.sub('_.+','_'+LocVar.get(re.search('_(.+)',x[0]).group(1)),x[0])+'" stoichiometry="'+str(x[1])+'" constant="true"/>\n' for x in P]# ID + Stoichometry
# ID should be MAM
# 	RxnProduct =  ['<speciesReference species="M_'+x[2]+'" stoichiometry="'+str(x[0])+'"></speciesReference>' for x in Rxn.Product] # ID + Stoichometry
		ProdRef = ''
		z = 0
		#print(SubsRef,SubsRef[10:].strip()) 
		#print(ProdRef,ProdRef[10:].strip())     
		while z < len(RxnProduct):
			ProdRef = ProdRef + RxnProduct[z]
			z = z+1
		#print(303, Rxn.GPR )        
		RxnGPRIdent = Rxn.GPR 
		if RxnGPRIdent[0].replace("[","").replace("]",""):    
#			RxnGPR = ['<modifierSpeciesReference species="E_'+x+'"/>' for x in re.findall('\[([A-Za-z0-9*]+)\]', Rxn.GPR[1])]
			RxnGPR = ['<modifierSpeciesReference species="'+x+'"/>' for x in re.findall('([A-Za-z0-9\-]+)', Rxn.GPR[1].replace("and","").replace("or",""))]
			#print(RxnGPR)    
			ModRef = ''
			z = 0
			while z < len(RxnGPR):
				ModRef = ModRef + RxnGPR[z]
				#print(RxnGPR)
				z = z+1 
		#print(ModRef)  
		try:        
			RxnGpr = Rxn.GPR2()[CL]
			#print(5, RxnGpr)            
			r30 =  RxnGpr       
			g2 =list(set(RxnGpr.split(' ')[1::2]))
			#print(g2) 
			if len(RxnGpr.split()) > 1 and g2:          
				g = '<fbc:geneProductAssociation>\n          <fbc:'+''.join([i for i in g2])+'>\n'+'\n'.join(['            <fbc:geneProductRef fbc:geneProduct="'+i+'"/>' for i in r30.split(' ')[::2]])+'\n          </fbc:'+''.join([i for i in g2])+'>\n        </fbc:geneProductAssociation>'
				#g = g.replace('<fbc:>','').replace('</fbc:>','')            
				#print(g)
			else: g='''
        <fbc:geneProductAssociation>
            <fbc:geneProductRef fbc:geneProduct="'''+RxnGpr+'''"/>
        </fbc:geneProductAssociation>    
        '''
			if not RxnGpr:
				#g = re.sub('<.?fbc:>','',g)
				#g = re.sub('  <fbc:geneProductRef','<fbc:geneProductRef',g) 
				g=''                
				#print(g)                
		except Exception:
			g=''    
		#print(333)
		no = mar_number(RxnID)        
		ModelRxn ='''      <reaction metaid="'''+no+'''" sboTerm="SBO:0000176" id="'''+no+'''" reversible="'''+RxnRev+'''" fast="false" fbc:lowerFluxBound="'''+lb+'''" fbc:upperFluxBound="FB3N1000">
        <notes>
          <body xmlns="http://www.w3.org/1999/xhtml">
            <p>Confidence Level: 0</p>
          </body>
        </notes>
        <annotation>
         <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:vCard4="http://www.w3.org/2006/vcard/ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
            <rdf:Description rdf:about="#'''+no+'''">
              <bqbiol:is>
                <rdf:Bag>
                  <rdf:li rdf:resource="https://identifiers.org/kegg.reaction/'''+str(RxnID.split("_")[0])+'''"/>
                  '''+RxnEC+'''
                </rdf:Bag>
              </bqbiol:is>
            </rdf:Description>
          </rdf:RDF>
        </annotation>
        <listOfReactants>
          '''+SubsRef[10:].strip()+'''
        </listOfReactants>
        <listOfProducts>
          '''+ProdRef[10:].strip()+'''
        </listOfProducts>
        '''+g+'''        
      </reaction>'''
		ModelRxn = re.sub(r'\n\s*\n','\n',ModelRxn)        
		return ModelRxn
	except Exception as error:
		print(error)
		return ''       
    
def objective():
	objective_reaction = '''    <fbc:listOfObjectives fbc:activeObjective="obj">
      <fbc:objective fbc:id="obj" fbc:type="maximize">
        <fbc:listOfFluxObjectives>
          <fbc:fluxObjective fbc:reaction="" fbc:coefficient="1"/>
        </fbc:listOfFluxObjectives>
      </fbc:objective>
    </fbc:listOfObjectives>'''
	return objective_reaction  

############### Modifier ############### 
def DefMod(Gene):  
	if Gene.Name():
		GeneID = Gene.Name().replace("[","").replace("]","").replace('-','')
		ID = '''<rdf:li rdf:resource="https://identifiers.org/hgnc.symbol/'''+GeneID+'''"/>'''  
	else:
		ID = ''
	if Gene.Entrez():
		GeneEntrez = Gene.Entrez()
		GeneEntrez = '''<rdf:li rdf:resource="https://identifiers.org/ncbigene/'''+str(GeneEntrez)+'''"/>\n'''
	else:
		GeneEntrez = ''
	if Gene.Ensg():
		GeneEnsbl = Gene.Ensg()
		try:        
			z = str(urllib.request.urlopen('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+GeneEnsbl).read())
			y = ['<rdf:li rdf:resource="https://identifiers.org/ensembl/'+str(x)+'"/>' for x in list(set(re.findall('ENST[0-9]+',z)))]
			y2 = ['<rdf:li rdf:resource="https://identifiers.org/ensembl/'+str(x)+'"/>' for x in list(set(re.findall('ENSP[0-9]+',z)))]      
		except Exception: # URLError: <urlopen error [Errno 60] Operation timed out>
			y, y2 = '',''            
	else:
		GeneEnsbl = Gene.Name()
		y = ''
		y2 = ''        
	ModelMod = '''      <fbc:geneProduct metaid="'''+GeneEnsbl+'''" sboTerm="SBO:0000243" fbc:id="'''+GeneEnsbl+'''" fbc:label="'''+GeneEnsbl+'''">
        <annotation>
          <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:vCard4="http://www.w3.org/2006/vcard/ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
            <rdf:Description rdf:about="#'''+GeneEnsbl+'''">
              <bqbiol:is>
                <rdf:Bag>
                  <rdf:li rdf:resource="https://identifiers.org/ensembl/'''+GeneEnsbl+'''"/>
                  '''+'\n                  '.join([str(x) for x in y])+'''
                  '''+'\n                  '.join([str(x) for x in y2])+'''
                  '''+GeneEntrez+'''
                  '''+ID+'''
                </rdf:Bag>
              </bqbiol:is>
            </rdf:Description>
          </rdf:RDF>
        </annotation>
      </fbc:geneProduct>'''
	ModelMod = re.sub(r'\n\s*\n','\n',ModelMod)         
	return ModelMod


def DefGrp(PathNameRxn): 
	listOfGroups=""
	for PathName in PathNameRxn:       
		groups=""
		group= '''      <groups:group sboTerm="SBO:0000633" groups:id="group" groups:name="'''+PathName+'''" groups:kind="partonomy">
        <groups:listOfMembers>\n'''      
		for RxnID in PathNameRxn[PathName].split():
			idRef = '''          <groups:member groups:idRef="'''+str(RxnID.split("_")[0])+'''"/>\n'''
			groups += idRef            
		listOfMembers= '''        </groups:listOfMembers>
      </groups:group>\n'''
		Group = group+groups+listOfMembers
		if not re.search('''<groups:listOfMembers>
        </groups:listOfMembers>''',Group,re.DOTALL):               
			listOfGroups += Group
		else:print(PathName)             
	#print(listOfGroups)        
	return listOfGroups
