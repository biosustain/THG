#!/usr/bin/python
# -*- coding: utf-8 -*-

from functions_generate_database import *
import re


class pathway(object):
	def __init__(self,url,time,ID,urlReferer,PathName):
		self.pagina = getHtml(url,time,urlReferer)
		if type(self.pagina) == bytes: self.pagina = self.pagina.decode('utf-8')
		self.link = getLinkPath(self.pagina)
		self.ID = ID
		self.PathName = PathName
	def ID(self): #KEGG
		try:		
			return self.ID
		except Exception:
			return ''
	def PathName(self): #KEGG
		try:		
			return self.PathName
		except Exception:
			return ''
	def Compounds(self): #KEGG
		try:		
			return [x for x in self.link[0]]
		except Exception:
			return ''
	def Reactions(self): #KEGG
		try:		
			return [x for x in self.link[1]]
		except Exception:
			return ''     


class reaction(object):
	def __init__(self,url,time,ID,path,termdyn,*newparam):
		self.pagina = getHtml(url,time)    
		self.link = getReacParam(self.pagina,time)
		self.ID = ID
		self.path = path
		self.termdyn = termdyn
		self.newparam = newparam  
	def ID(self): #Patway-KEGG
		try:		
			return self.ID
		except Exception:
			return ''      
	def Name(self): #KEGG
		try:	
			if self.link[3]:
				N = self.link[3]
			else:
				N = self.ID
			return N
		except Exception:
			return ''
	def EC(self): #KEGG
		try:		
#			return self.link[4][0][1]
			return [x[1] for x in self.link[4]]
		except Exception:
			return ''
	def GPR(self): #MetaCyc
		try:           
			return self.newparam[0][0] , self.newparam[0][1]
		except Exception:
			return ''
#	def GPR2(self): #MetaCyc
#		try:		          
#			return self.link[8]
#		except Exception:
#			return ''        
	def Termodyn(self): #Patway-KEGG
		try:				
			return self.termdyn
		except Exception:
			return ''
	def Substrate(self): #KEGG     
		try:		           
			S = [[x[0]]+[x[1]]+[x[2]] for x in self.link[1]]           
			return  S #link+ID+stoichometry
		except Exception:
			return ''
	def SetSubstrate(self, substrate): #KEGG'        
		self.link[1] = [[x[0]]+[x[1]]+[x[2]] for x in substrate]
		self.link[0] = self.link[1] + self.link[2]         
		return self.link[1]    
	def Product(self): #KEGG
		try:		
			P = [[x[0]]+[x[1]]+[x[2]] for x in self.link[2]]               
			return  P #link+ID+stoichometry
		except Exception:
			return ''
	def SetProduct(self, product): #KEGG         
		self.link[2] = [[x[0]]+[x[1]]+[x[2]] for x in product]
		self.link[0] = self.link[1] + self.link[2]         
		return self.link[2]     
	def Pathway(self): #Patway-KEGG
		try:		
			return self.path
		except Exception:
			return ''
	def Subcel(self): #Uniprot
		try:					
			return  self.newparam[0][2] , self.newparam[0][3]
		except Exception:
			return ''
	def Equivalent(self): #KEGG
		try:					
			return  self.link[5]
		except Exception:
			return ''
	def GTest(self): #KEGG
		try:					
			return  self.link[6]
		except Exception:
			return ''
	def CTest(self): #KEGG
		try:					
			return  self.link[7]
		except Exception:
			return ''
	def MBTest(self): #KEGG
		try:					
			return  self.ID[0]
		except Exception:
			return ''


class gpr(object):
	def __init__(self,ec,time):
		self.ec = ec
		session = setup_biocyc_session()
		self.GPRPAss = getGPR(self.ec,session)  
		self.Subcell = getLocation(self.GPRPAss[3],self.GPRPAss[1],self.GPRPAss[2],1, 'files/bb.pickle',session)  
	def EC(self): #Patway-KEGG
		try:		
			return self.ec
		except Exception:
			return ''
	def GprSubcell(self): #MetaCyc 
		try:		           
			return self.GPRPAss[3] , self.GPRPAss[4] , self.Subcell[0] , self.Subcell[1] # General gene rule with stoichometry + General gene rule without stoichometry + Subcell-specific rule with stoichometry + Subcell-specific rule without stoichometry 
		except Exception:
			return ''


class gene(object):
	def __init__(self,gene,db):
		self.gene = gene
		self.db = db
	def Name(self): #Patway-KEGG
		try:		
			return self.gene
		except Exception:
			return ''
	def Ensg(self): #MetaCyc 
		try:		
			iiii2 =re.search('gene=([A-Z0-9]+)',str(urllib.request.urlopen('https://www.genome.jp/dbget-bin/www_bget?hsa+'+self.gene).read()))#.group(1) 
			if not iiii2:
				iiii2 = re.search('(ENSG[0-9]+)',str(urllib.request.urlopen('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+self.gene).read()))                     
			if iiii2:
				EnsGene = iiii2.group(1)
			return EnsGene
		except Exception:
			return ''
	def Entrez(self): #MetaCyc 
		try:		
			EntrezGene = re.findall(self.gene+"_HUMAN[\S\s].*?\n", open(self.db).read())[0].split('\t')[4]
			return EntrezGene
		except Exception:
			return ''
	def Uniprot(self): #Uniprot 
		try:		
			UniProtGene = re.findall(self.gene+"_HUMAN[\S\s].*?\n", open(self.db).read())[0].split('\t')[0]
			return UniProtGene
		except Exception:
			return ''


class compound(object):
	def __init__(self,url,ident,time,EF,specialCompounds,*newparam):
		self.ident = ident
		self.pagina = str(urllib.request.urlopen('https://www.genome.jp/entry/'+self.ident).read()) #str(getHtml(url,20))
		self.atributes = getCompParam(self.pagina,self.ident,time,EF,specialCompounds, RxnID=None)    
		self.newparam = newparam
	def ID1(self):
		try:		           
			return self.atributes[0][0][0]
#			return getID(url)
		except Exception:
			return ''
	def ID2(self):
		try:		
			return self.atributes[0][0][1]
#			return getID(url)
		except Exception:
			return ''
	def AssRxn1(self):
		try:		
			return self.atributes[2][0]
		except Exception:
			return ''
	def AssRxn2(self):
		try:		
			return  [x[0] for x in self.atributes[2][0]]
		except Exception:
			return ''
	def AssRxn3(self):
		try:		
			return [x[1] for x in self.atributes[2][0]] 	
		except Exception:
			return ''
	def Formula1(self):
		try:	          
			return self.atributes[1][0][0]
		except Exception:
			return ''
	def Formula2(self):
		try:	            
			return self.atributes[1][0][1]
		except Exception:
			return ''
	def Formula3(self):
		try:	            
			return self.atributes[1][0][0]
		except Exception:
			return ''
	def Atom1(self):
		try:
			if self.atributes[0][0][0][0] == 'C':
				composition = atom(self.atributes[1][0][0])
			elif self.atributes[0][0][0][0] == 'G':
				composition = glycan(self.atributes[1][0][0])
			return composition
		except Exception:
			return ''
	def Atom2(self):
		try:
			if self.atributes[0][0][1][0] == 'C':
				composition = atom(self.atributes[1][0][1])
			elif self.atributes[0][0][1][0] == 'G':
				composition = glycan(self.atributes[1][0][1])
			return composition
		except Exception:
			return ''
	def Atom3(self):
		try:
			if self.atributes[0][0][0][0] == 'C':               
				composition = atom(self.atributes[1][0][0])          
			elif self.atributes[0][0][0][0] == 'G':  
				composition = glycan(self.atributes[1][0][0])             
			return composition
		except Exception:
			return ''
	def Name(self):
		try:		
			return self.atributes[3][0]
		except Exception:
			return ''
	def Subcel(self):
		try:		
			return self.newparam
		except Exception:
			return ''
	def PubChem(self):
		try:		           
			return self.atributes[4][0]
		except Exception:
			return ''
	def CheBI(self):
		try:		
			return self.atributes[5][0]
		except Exception:
			return ''
	def LIPIDMAPS(self):
		try:		
			return self.atributes[6][0]
		except Exception:
			return ''
	def LipidBank(self):
		try:		
			return self.atributes[7][0]
		except Exception:
			return ''
	def GlyDB(self):
		try:		
			return self.atributes[8][0]
		except Exception:
			return ''
	def JCGGDB(self):
		try:		
			return self.atributes[9][0]
		except Exception:
			return ''
	def charge(self):
		try:		           
			return self.atributes[10][:]
		except Exception:
			return ''
	def Formula4(self):
		try:		           
			return self.atributes[11][:]
		except Exception:
			return ''
	def inchikey(self):
		try:		           
			return self.atributes[12][:]
		except Exception:
			return ''
	def inchi(self):
		try:		           
			return self.atributes[13][:]
		except Exception:
			return ''
	def CID(self):
		try:          
			return self.atributes[14][:]
		except Exception as e:          
			return ''  
