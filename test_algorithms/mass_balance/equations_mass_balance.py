###########
# This piece of code is fully functional and ready to be integrated into MB_REM function to choose the solution with the smallest steochiometric coefficients. The code is able to use the output of solve function (line 30)
# In addition a new function has been created "maximumGCD" (line 58) that has to be added to the code in function.py
from cobra.io import read_sbml_model
from lib2to3.pgen2.token import GREATER
from cmath import isnan
import sys
import numpy as np
import scipy
#import sys
import re
import json
import csv
from collections import defaultdict
#from scipy.optimize import fsolve
#from scipy.optimize import nnls 
#from scipy.sparse.linalg import lsqr
from sympy.solvers import solve
#from sympy import Symbol
#import matplotlib.pyplot as plt
#from fractions import gcd
import numpy.matlib
#from sympy import *
import sympy
from numpy import size#,matrix,zeros
#from numpy.linalg import det
#import scipy.optimize as optimize
from scipy.optimize import linprog
from itertools import chain, zip_longest 
import urllib.request
from functools import reduce
import traceback
#from math import gcd
import math

### Required functions called by MB functions:
### Compound: Determine the atomic composition of a compound
def atom10(Formula):
	try:          
		if 'C' in Formula:
			C = re.findall(r'C([0-9]+)', Formula)
			if not C:
				C = ['1']
		else:
			C = ['0']
		if 'H' in Formula:
			H = re.findall(r'H([0-9]+)', Formula)
			if not H:
				H = ['1']
		else:
			H = ['0']
		if 'O' in Formula:
			O = re.findall(r'O([0-9]+)', Formula)
			if not O:
				O = ['1']
		else:
			O = ['0']
		if 'N' in Formula:
			N = re.findall(r'N([0-9]+)', Formula)
			if not N:
				N = ['1']
		else:
			N = ['0']
		if 'P' in Formula:
			P = re.findall(r'P([0-9]+)', Formula)
			if not P:
				P = ['1']
		else:
			P = ['0']
		if 'S' in Formula:
			S = re.findall(r'S([0-9]+)', Formula)
			if not S:
				S = ['1']
		else:
			S = ['0']
		if 'K' in Formula:
			K = re.findall(r'K([0-9]+)', Formula)
			if not K:
				K = ['1']
		else:	
			K = ['0']
		if 'Ca' in Formula:
			Ca = re.findall(r'Ca([0-9]+)', Formula)
			if not Ca:
				Ca = ['1']
		else:	
			Ca = ['0']
		if 'Na' in Formula:
			Na = re.findall(r'Na([0-9]+)', Formula)
			if not Na:
				Na = ['1']
		else:	
			Na = ['0']
		if 'Fe' in Formula:
			Fe = re.findall(r'Fe([0-9]+)', Formula)
			if not Fe:
				Fe = ['1']
		else:	
			Fe = ['0']
		if 'F' in Formula:
			F = re.findall(r'F([0-9]+)', Formula)
			if not Fe:
				F = ['1']
		else:	
			F = ['0']            
		if 'X' in Formula:
			X = re.findall(r'X([0-9]+)', Formula)
			if not X:
				X = ['1']
		else:	
			X = ['0']
		if 'R' in Formula:
			R = re.findall(r'R([0-9]+)', Formula)
			if not R:
				R = ['1']
		else:	
			R = ['0']			      
		return [int(C[0]), int(H[0]), int(O[0]), int(N[0]), int(P[0]), int(S[0]), int(K[0]), int(Ca[0]), int(Na[0]), int(Fe[0]), int(X[0]), int(R[0])]
	except Exception as e:        
		return C, H, O, N, P, S, K, Ca, Na, Fe, F, X, R

    
### Required functions called by MB functions:
### Compound: Determine the atomic composition of a compound
def atom100(Formula):
	try:          
		if 'C' in Formula:
			C = re.findall(r'(C[0-9]+)', Formula)
			if not C and re.findall('C(?![a-z]+)', Formula):
				C = ['C1']
		else:
			C = [] 
		if 'H' in Formula:
			H = re.findall(r'(H[0-9]+)', Formula)
			if not H and re.findall('H(?![a-z]+)', Formula):
				H = ['H1']
		else:
			H = [] 
		if 'O' in Formula:
			O = re.findall(r'(O[0-9]+)', Formula)
			if not O and re.findall('O(?![a-z]+)', Formula):
				O = ['O1']
		else:
			O = [] 
		if 'N' in Formula:
			N = re.findall(r'(N[0-9]+)', Formula)             
			if not N and re.findall('N(?![a-z]+)', Formula):
				N = ['N1']
		else:
			N = [] 
		if 'P' in Formula:
			P = re.findall(r'(P[0-9]+)', Formula)
			if not P and re.findall('P(?![a-z]+)', Formula):
				P = ['P1']
		else:
			P = []
		if 'Se' in Formula:
			Se = re.findall(r'(Se[0-9]+)', Formula)
			if not Se:
				Se = ['Se1']
		else:	
			Se = []             
		if 'S' in Formula:
			S = re.findall(r'(S[0-9]+)', Formula)           
			if not S and re.findall('S(?![a-z])', Formula):            
				S = ['S1']
		else:
			S = [] 
		if 'K' in Formula:
			K = re.findall(r'(K[0-9]+)', Formula)
			if not K and re.findall('K(?![a-z]+)', Formula):
				K = ['K1']
		else:	
			K = [] 
		if 'Ca' in Formula:
			Ca = re.findall(r'(Ca[0-9]+)', Formula)
			if not Ca:
				Ca = ['Ca1']
		else:	
			Ca = [] 
		if 'Na' in Formula:
			Na = re.findall(r'(Na[0-9]+)', Formula)
			if not Na:
				Na = ['Na1']
		else:	
			Na = [] 
		if 'Fe' in Formula:
			Fe = re.findall(r'(Fe[0-9]+)', Formula)
			if not Fe:
				Fe = ['Fe1']
		else:	
			Fe = [] 
		if 'F' in Formula:
			Fe = re.findall(r'F([0-9]+)', Formula)
			if not F and not re.findall('F(?![a-z]+)', Formula):
				F = ['1']
		else:	
			F =  []           
		if 'X' in Formula:
			X = re.findall(r'(X[0-9]+)', Formula)
			if not X and re.findall('X(?![a-z]+)', Formula):
				X = ['X1']
		else:	
			X = [] 
		if 'R' in Formula:
			R = re.findall(r'(R[0-9]+)', Formula)
			if not R and re.findall('R(?![a-z]+)', Formula):
				R = ['R1']
		else:	
			R = []    
		if 'M' in Formula:
			M = re.findall(r'(M[0-9]+)', Formula)
			if not M and re.findall('M(?![a-z]+)', Formula):
				M = ['M1']
		else:	
			M = []
		if 'I' in Formula:
			I = re.findall(r'(I[0-9]+)', Formula)
			if not I and re.findall('I(?![a-z]+)', Formula):
				I = ['I1']
		else:	
			I = []             
		return ''.join(str(i) for i in C+H+O+N+P+S+K+Ca+Na+Fe+F+X+M+Se+I+R)
	except Exception as e:        
		return Formula 
    

### Transform to np array
def inarray(A,B):
    try:
        A=np.array(A)
        B=np.array(B)        
		# remove indices where A is null
        m = A!=0
		# ensure B values are null where A value are null
        assert np.allclose(B[~m], 0)
		# compute B/A
        out = B[m]/A[m]
		# ensure all values are (almost) equal
        assert np.allclose(out-out[0], 0)
		# print result
        if (out[0]).is_integer() and out[0]>0:
            return int(out[0])
    except Exception as e:
        print(e)
        return ''

## Transform equation to matrix 
def eq2mat(eq):  
	#Transform Equation to Matrix
	Subs = [x.split(' + ')  for x in eq.split(' -> ')][0]
	Prod = [x.split(' + ')  for x in eq.split(' -> ')][1]
	Met = Subs+Prod
	a = [re.findall('[A-Z]',x) for x in Met]
	AtomList = []
	for j in a:
		for x in j:
			if not x in AtomList:
				AtomList = AtomList + [x]
	Matrix = numpy.zeros(shape=(len(AtomList),len(Met)))
	i = 0
	while i < len(AtomList): #loop for atoms (row)	
		j = 0
		while j < len(Met): #loop for compounds (column)
			AtomJth = re.findall(AtomList[i]+'([0-9]+)',Met[j])
			if not AtomJth:
				AtomJth = re.findall(AtomList[i],Met[j])
				if AtomJth:
					AtomJth = [1]
				else:
					AtomJth = [0]
#			if Met[j] in Prod: 
#				Matrix[i,j] = str(int(AtomJth[0])*(-1))
#			else:
#				Matrix[i,j] = AtomJth[0]
			Matrix[i,j] = AtomJth[0]
			j += 1
		i += 1	       
	return Matrix

## Define null space 
def nullity(Matrix):
#	print(Matrix)    
	Matrix=Matrix.astype(float)
	Matrix2 =Matrix.tolist(  ) 
	#Calculate the nullity or dimentionMallity
	rank_m = np.linalg.matrix_rank(Matrix2, tol=None)  
	nullity = len(Matrix2[0])-rank_m
	if nullity != 0 or len(Matrix2[0]) != len(Matrix2):
		#Calculate Linearly independent Matrix
		rank_ref = 1
		R_independent = Matrix[0]
		R2 = Matrix[0]
		i=1
		while i < len(Matrix2): # antes Matrix2[0]
			R2 = np.vstack([R2, Matrix[i,:]])
			rank_ith = np.linalg.matrix_rank(R2, tol=None)
			if rank_ith <= rank_m and rank_ith != rank_ref:
				rank_ref = rank_ith
				R_independent = np.vstack([R_independent, Matrix[i,:]])
			i += 1
		#Extended matrix
		ExpMatrix = [[0 for col in range(len(Matrix2[0]))] for row in range(nullity)]
		if ExpMatrix:
			c=0
			while c < np.shape(ExpMatrix)[0]:
				g = 0
				while g <= c:
					ExpMatrix[c][len(Matrix2[0])-(g+1)]=1
					g = g+1
				c = c+1
			ExpMatrix2 = np.vstack([R_independent, ExpMatrix])
		else:
			ExpMatrix2 = R_independent
	else:
		ExpMatrix2 = Matrix
		R_independent = Matrix        
	return ExpMatrix2 , R_independent

## Inverse Matrix
def inv(A):
	#Cofactor Matrix
	MC = numpy.zeros(shape=(len(A),len(A)))
	r = 1
	for i in range(size(A[0])):
		for j in range(size(A[1])):
			cof = scipy.delete(scipy.delete(A, i, 0), j, 1)
			MC[i,j] = numpy.linalg.det(cof)*r
			if abs(MC[i,j]) < 0.001:
				MC[i,j] = 0.
			r = -1*r
	#Adjunt Matrix
	MADJ=MC.transpose() # Matriz adjunta 
	#Determinant
	determ = numpy.linalg.det(A)
	#Inverse Matrix
	IM = 1/determ*(MADJ) # inverse_matrix = 1/det(A)*adjunt_matrix where adjunt_matrix = transposed_cofactor_matrix
	#np.linalg.inv(A)    
	return IM

## New: Calculates the greatest common denominator
def maximumGCD(A, I, K):
    # Calculates the overall gcd between all the variables of an undetermined system for a given value of the idependent term
    # A= Array, I=Independent term, K= Value of the independent term     
    N=len(A)
    # Initialize a dp table of size N*2
    dp = [0 for x in range(N*N)]
    for j in range(0,N):
        Jth = int(eval(str(A[j]).replace("'","").replace(I, "K")))
        # Traverse the array A[] over indices [1, N - 1]
        for i in range(j+1, N):
            # Store the previous state results
            Ith = int(eval(str(A[i]).replace("'","").replace(I, "K")))
            # Store maximum GCD result
            # for each current state
            dp[j*N+i] = math.gcd(Jth, Ith)        
    dp=dp+np.where(np.array(dp), np.array(dp)==0, max(dp))
    mx=min(dp)	             
    # Return the result
    return mx 

## Identify differences between reactions
def RxnCompare(eq1,eq2):
	SpeciesInEq1 = ([x for x in [x.replace(" ","") for x in eq1.split('->')][0].split('+')],[x for x in [x.replace(" ","") for x in eq1.split('->')][1].split('+')])
	SpeciesInEq2 = ([x for x in [x.replace(" ","") for x in eq2.split('->')][0].split('+')],[x for x in [x.replace(" ","") for x in eq2.split('->')][1].split('+')])
	NewSpeciesInSubstrates = [x for x in SpeciesInEq1[0] if x not in SpeciesInEq2[0]]
	NewSpeciesInProducts = [x for x in SpeciesInEq1[1] if x not in SpeciesInEq2[1]]
	IndexS = [-1 for x in SpeciesInEq1[0] if x not in SpeciesInEq2[0]]
	IndexP = [1 for x in SpeciesInEq1[1] if x not in SpeciesInEq2[1]]
	NewSpecies = NewSpeciesInSubstrates+NewSpeciesInProducts
	NewIndex = IndexS+IndexP
	return NewSpecies, NewIndex

### Glycans: Transform a glycan reaction to be MB 
def Reformulation(eq):       
	from Function import glycan
	CompLs = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
	#define a dictionary
	GLIdent = []
	GlList = {}
	p = 0
	while p < len(eq.split('->')):
		k = 0
		while k < len(eq.split('->')[p].split('+')):
			Gl = glycan(eq.split('->')[p].split('+')[k])
			g = 0
			while g < len(Gl):
				if not Gl[g][0][0] in GLIdent:
					GLIdent = GLIdent + [Gl[g][0][0]]
					GlList[Gl[g][0][0]] = CompLs.pop(0)				
				g = g + 1 			
			k = k +1 
		p = p + 1
	#Transform Substrates notation
	Sb = ''
	p = 0
	while p < 1:
		k = 0
		while k < len(eq.split('->')[p].split('+')):
			GlF = glycan(eq.split('->')[p].split('+')[k])
			g = 0
			c = ''
			while g < len(GlF):
				a = GlF[g][0][0].replace(GlF[g][0][0],GlList[GlF[g][0][0]])
				b = GlF[g][0][1]
				c = c+str(a)+str(b)
				g = g + 1
			Sb = Sb+' + '+c
			k = k + 1
		p = p + 1				
	#Transform Products notation
	Pr = ''
	p = 1
	while p < 2:
		k = 0
		while k < len(eq.split('->')[p].split('+')):
			GlF = glycan(eq.split('->')[p].split('+')[k])
			g = 0
			c = ''
			while g < len(GlF):
				a = GlF[g][0][0].replace(GlF[g][0][0],GlList[GlF[g][0][0]])
				b = GlF[g][0][1]
				c = c+str(a)+str(b)
				g = g + 1
			Pr = Pr+' + '+c
			k = k + 1
		p = p + 1
	# Transform reaction notation
	SP = Sb[3:]+' -> '+Pr[3:]    
	return SP

### Extract the inital parameters from Reaction atributes
def WrapRxnSubsProdParam(Reaction,MetList,MetEquiv): 
	Substrate = Reaction.Substrate()
	Product = Reaction.Product()
	DS = {} # Initial empty diccionary Substrates
	DP = {} # Initial empty diccionary Products
	# Generate the indexes in the dictionary
	if str(Reaction.GTest()) == '1':
		for s in range(0,len(Substrate)):
			if Substrate[s][2] in MetEquiv:
				DS[MetList[MetEquiv[Substrate[s][2]]].Formula2()] = [1 , MetList[MetEquiv[Substrate[s][2]]].ID1()]  #Formula1
			else:
				DS[MetList[Substrate[s][2]].Formula2()] = [1 , MetList[Substrate[s][2]].ID1()] #Formula1
		for p in range(0,len(Product)):
			if Product[p][2] in MetEquiv:
				DP[MetList[MetEquiv[Product[p][2]]].Formula2()] = [1 , MetList[MetEquiv[Product[p][2]]].ID1()]  #Formula1
			else:
				DP[MetList[Product[p][2]].Formula2()] = [1 , MetList[Product[p][2]].ID1()] #Formula1
	else:
		for s in range(0,len(Substrate)):
			if Substrate[s][2] in MetEquiv:
				DS[MetList[MetEquiv[Substrate[s][2]]].Formula1()] = [1 , MetList[MetEquiv[Substrate[s][2]]].ID1()]  #Formula1
			else:
				DS[MetList[Substrate[s][2]].Formula1()] = [1 , MetList[Substrate[s][2]].ID1()] #Formula1
		for p in range(0,len(Product)):
			if Product[p][2] in MetEquiv:
				DP[MetList[MetEquiv[Product[p][2]]].Formula1()] = [1 , MetList[MetEquiv[Product[p][2]]].ID1()]  #Formula1
			else:
				DP[MetList[Product[p][2]].Formula1()] = [1 , MetList[Product[p][2]].ID1()] #Formula1
	return DS, DP

### Incorporates the results of the MB process to reaction
def UnwrapRxnSubsProdParam(Reaction, LibIni, IthRxnMB):
	#Substrates
	SubsLibEnd = dict(x for x in list(zip(IthRxnMB[6],IthRxnMB[0])))
	for x in SubsLibEnd: 
		if x in LibIni[0]:
			LibIni[0][x][0] = SubsLibEnd[x] #Change a default St Coef (1) by the calculated one
		else: # Metabolite in the new list but not in the original =  addedd metabolite 
			if x == 'Fe':
				LibIni[0].update({x:[SubsLibEnd[x],'C00023']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'Na':
				LibIni[0].update({x:[SubsLibEnd[x],'C01330']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'Ca':
				LibIni[0].update({x:[SubsLibEnd[x],'C00076']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'K':
				LibIni[0].update({x:[SubsLibEnd[x],'C00238']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'F':
				LibIni[0].update({x:[SubsLibEnd[x],'C00023']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'R':
				LibIni[0].update({x:[SubsLibEnd[x],'C00000']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'X':
				LibIni[0].update({x:[SubsLibEnd[x],'C0000X']}) # This line has to be replaced by a function that gives the real ID
	ProdLibEnd = dict(x for x in list(zip(IthRxnMB[7],IthRxnMB[1])))
	for x in ProdLibEnd: 
		if x in LibIni[1]:
			LibIni[1][x][0] = ProdLibEnd[x] #Change a default St Coef (1) by the calculated one
		else: # Metabolite in the new list but not in the original =  addedd metabolite 
			if x == 'Fe':
				LibIni[0].update({x:[SubsLibEnd[x],'C00023']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'Na':
				LibIni[0].update({x:[SubsLibEnd[x],'C01330']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'Ca':
				LibIni[0].update({x:[SubsLibEnd[x],'C00076']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'K':
				LibIni[0].update({x:[SubsLibEnd[x],'C00238']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'F':
				LibIni[0].update({x:[SubsLibEnd[x],'C00023']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'R':
				LibIni[0].update({x:[SubsLibEnd[x],'C00000']}) # This line has to be replaced by a function that gives the real ID
			elif x == 'X':
				LibIni[0].update({x:[SubsLibEnd[x],'C0000X']}) # This line has to be replaced by a function that gives the real ID
	return LibIni

### Called directly by RxnBalance functions:
## 1. Extracts the reaction formula in the form of string from the parameters of an Object of Class Reaction. Only to create DB
def RxnParam2Eq(Reaction,MetList,MetEquiv):  
	Substrate = Reaction.Substrate()
	Product = Reaction.Product()
	if str(Reaction.GTest()) == '1':
		s = 0
		S = ''
		while s < len(Substrate):  
			SStch = str(Substrate[s][0])
			if Substrate[s][2] in MetEquiv:                   
				S = S+' + '+SStch+' '+MetList[MetEquiv[Substrate[s][2]]].Formula2() #Formula2
			else:
				S = S+' + '+SStch+' '+MetList[Substrate[s][2]].Formula2() #Formula2
			s = s + 1
		S = S[3:].replace(" +  + "," + ")
		p = 0
		P = ''
		while p < len(Product):
			PStch = str(Product[p][0])
			if Product[p][2] in MetEquiv:
				P = P+' + '+PStch+' '+MetList[MetEquiv[Product[p][2]]].Formula2() #Formula2 
			else:
				P = P+' + '+PStch+' '+MetList[Product[p][2]].Formula2() #Formula2
			p = p + 1
		P = P[3:].replace(" +  + "," + ")
		eq0 = S+" -> "+P          
		eq = Reformulation(eq0)
	if str(Reaction.GTest()) != '1':
		s = 0
		S = ''
		while s < len(Substrate):
			SStch = str(Substrate[s][0])
			if Substrate[s][2] in MetEquiv:
				S = S+' + '+SStch+' '+MetList[MetEquiv[Substrate[s][2]]].Formula1()
			else:
				S = S+' + '+SStch+' '+MetList[Substrate[s][2]].Formula1()
			s = s + 1
		S = S[3:].replace(" +  + "," + ")
		p = 0
		P = ''
		while p < len(Product):
			PStch = str(Product[p][0])
			if Product[p][2] in MetEquiv:
				P = P+' + '+PStch+' '+MetList[MetEquiv[Product[p][2]]].Formula1() 
			else:
				P = P+' + '+PStch+' '+MetList[Product[p][2]].Formula1() 
			p = p + 1
		P = P[3:].replace(" +  + "," + ")           
		eq0 = S+" -> "+P
		eq = eq0  
	if re.findall('\)n', eq): # check if it is a general formula
		eq = re.sub('\)n','',re.sub('\(','',re.sub('\)n[A-Za-z0-9]+','',eq)))
	return eq

## 2. New and improved version: Adds missing atoms to equation to enable mass balance (Ca, Na, Fe, X, R, K)
def AddMissingAtom(eq): 
	Ls = list('abcdefghijklmnopqrstuvwxyz')
	#Mass Balance
	Ss,Es,i,ii=defaultdict(list),[],1,defaultdict(list)
	for p in eq.split('->'):
		for k in p.split('+'):           
			c = [Ls.pop(0), 1]         
			for e,m in re.findall('([A-Z][a-z]?)([0-9]*)',k):              
				m = 1 if m == '' else int(m)
				d = [c[0],c[1]*m*i]
				Ss[e][:0],Es[:0] = [d],[[e,d]]            
				ii[e].append(d[1]) 
		i=-1
	k = [x[0] for x in list(ii.items()) if x[0] in ['Fe', 'X', 'R', 'Na', 'K', 'Ca', 'F'] and len(x[1]) == 1]
	if k and len(re.split('\+|->',eq)) > 2 and len(re.split('->',eq)[1]) > 2:  
		for x in k:
			if ii[x][0] < 0:
				eq = re.sub(" ->", " + "+x+" ->", eq)
			else:
				eq = re.sub(" ->", " -> "+x+" +", eq) 
	return eq

## 3.1 New CountAtom.  Here we do not check if metabolites are already in the list of metabolites or are new and have to be added- Stoichimetric indexes are considered
def CountAtom(eq,AddH,H2O,RxnID):  
	eq=" " + eq # add an extra space at the beginning to avoid problems when defining StCoeff of the first metabolite using reduce func    
	Ls = list('abcdefghijklmnopqrstuvwxyz')
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
#	MassBalance = [round(sum(x)) for x in ii.values()]
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
		SubsSpc = [x.split(' ')[1] for x in [x.strip() for x in Species[0]]]
		ProdStc = [x.split(' ')[0] for x in [x.strip() for x in Species[1]]]
		ProdSpc = [x.split(' ')[1] for x in [x.strip() for x in Species[1]]]
		SubsStch = [float(x) for x in SubsStc]
		ProdStch = [float(x) for x in ProdStc]            
	else:
		SubsStch = ''
		ProdStch = ''   
	return SubsStch , ProdStch, eq

## 3.2 First option to MB 2. Once finished, this will replace the previous
def MB_Core(eq,AddH,H2O,RxnID):
	Ls = list('abcdefghijklmnopqrstuvwxyz')
	LS = len(eq.split('->')[0].split('+'))
	#Mass Balance
	Ss,Os,Es,a,i=defaultdict(list),Ls[:],[],1,1
	for p in eq.split('->'):
		for k in p.split('+'):
			c = [Ls.pop(0), 1]
			for e,m in re.findall('([A-Z][a-z]?)([0-9]*)',k):
				m = 1 if m == '' else int(m)
				a*= m
				d = [c[0],c[1]*m*i]
				Ss[e][:0],Es[:0] = [d],[[e,d]]
		i=-1
	Ys = dict((s,eval('sympy.Symbol("'+s+'")')) for s in Os if s not in Ls)
	Qs = [eval('+'.join('%d*%s'%(c[1],c[0]) for c in Ss[s]),{},Ys) for s in Ss]+[Ys['a']-a]
	k = solve(Qs,*Ys)
	BadKTest = [x for x in [str(k[x]<=0) for x in [x for x in k]] if "True" in x]
	if k:
		I = set(re.findall('[a-z]', str(list(k.values())))) # Independent variable
		TestI = len(I) # Test if more than 1 independent variable
	else:
		I = ''
		TestI = 2	
	if k and not BadKTest and TestI==0:
		N = [k[Ys[s]] for s in sorted(Ys)]
		g = N[0]
		for a1, a2 in zip(N[0::2],N[1::2]):
			print(g,a2)            
			g = math.gcd(int(g),int(a2))
		N = [x/g for x in N]
		SubsStch = N[0:LS]
		ProdStch = N[LS:len(N)]
	elif k and not BadKTest and TestI==1:
		I = str(I).replace("{","").replace("}","").replace("'","")
		A = list(k.values())+[I] # List of variables' values. Here the independent variable is in the last position
		a = str(list(filter(lambda p: not I in p, str(A).replace("[","").replace("]","").split(',')))).replace("'","").replace(" ","")
		V = 1 # Variation
		Kref = int(min(np.matrix(a).A[0]))
		Solution = np.array(np.matrix(eval(str(A).replace("'","").replace(I,str(Kref))))).ravel()
#		Solution = [numpy.round(x,1) for x in abs(Solution)/min(abs(Solution))]
		TestSol = sum([1 for x in Solution if x <=0]) # Check for possible negative values
		if TestSol == 0:
			Solution = [numpy.round(x,1) for x in Solution/min(abs(Solution))]
			LS = len(eq.split('->')[0].split('+'))
			SubsStch = Solution[0:LS]
			ProdStch = Solution[LS:len(Solution)]
		else:
			SubsStch = ''
			ProdStch = ''
	else: # The system has no solution or multiple solutions with more than one indep. variable
		SubsStch = ''
		ProdStch = ''
	return SubsStch, ProdStch

## 3.3 New MB_REM. This version works well: weaknesses: 1. improve the exploration os sos, ii. only works when only 1 indep var
def MB_REM(eq,AddH,H2O,RxnID):    
	M = eq2mat(eq)  # equaction to matrix  
	A = nullity(M)	# determine the null space 
	# determine the inverse matrix in order to get the solution
	Solution = inv(A[0])[:,len(A[0])-1] # the last column of the inverse matrix A
	# If not solution then reduce original Matrix to Row-echelon form
	BadSolTest_a = [1 if i == 0 else 0 for i in Solution.tolist()] # Any stoichimetric coefficient = 0?
	BadSolTest_b = [1 if not str(Solution.tolist()).isdigit() else 0 for i in Solution.tolist()]# any stoichimetric coeficient is dependent?
	BadSolTest = sum(BadSolTest_a+BadSolTest_b)
	if not Solution.all() or BadSolTest!=0:
		REM = np.array(sympy.Matrix(M).rref()[0])
		REM = REM[~np.all(REM == 0, axis=1)]
		A = nullity(REM)
		Solution = inv(A[0])[:,len(A[0])-1]	
	BadSolTest_a = [1 if i == 0 else 0 for i in Solution.tolist()] # Any stoichimetric coefficient = 0?
	BadSolTest_b = [1 if not str(Solution.tolist()).isdigit() else 0 for i in Solution.tolist()]# any stoichimetric coeficient is dependent?
	BadSolTest = sum(BadSolTest_a+BadSolTest_b)
	if not Solution.all() or BadSolTest!=0:
		## Define the system's undetermined solution: copy paste from MB_Core
		Ls = list('abcdefghijklmnopqrstuvwxyz')
		LS = len(eq.split('->')[0].split('+'))
		#Mass Balance
		Ss,Os,Es,a,i,ii=defaultdict(list),Ls[:],[],1,1,defaultdict(list)
		for p in eq.split('->'):
		    for k in p.split('+'):
		        c = [Ls.pop(0), 1]
		        for e,m in re.findall('([A-Z][a-z]?)([0-9]*)',k):              
		            m = 1 if m == '' else int(m)
		            a*= m
		            d = [c[0],c[1]*m*i]
		            Ss[e][:0],Es[:0] = [d],[[e,d]]            
		            ii[e].append(d[1])
		    i=-1    
		Ys = dict((s,eval('sympy.Symbol("'+s+'")')) for s in Os if s not in Ls)
		Qs = [eval('+'.join('%d*%s'%(c[1],c[0]) for c in Ss[s]),{},Ys) for s in Ss]+[Ys['a']-a]
		k = solve(Qs,*Ys)
		## Explore the space of soultions in order to find the solution with the overall smallest stoichimetric coefficients: new part
		if k:
			I = re.search('[a-z]', str(list(k.values()))) # Independent variable
			TestI = len(set(re.findall('[a-z]', str(list(k.values())))))# Test if more than 1 independent variable
		else:
			I = ''
			TestI = 0
		if I and TestI == 1: # If more than 1 indep variable, skip
			I = I.group()
			A = list(k.values())+[I] # List of variables' values. Here the independent variable is in the last position
			a = str(list(filter(lambda p: not I in p, str(A).replace("[","").replace("]","").split(',')))).replace("'","").replace(" ","")
			Kref = int(min(np.matrix(a).A[0]))# +V # Value of the independent term. In order to minimize the iteration task Kref has the value of the highest non-dependent 			variable and is -1 to enable a minimum iteration
			V = Kref/1000 # Variation
			if V < 1: V = 1
			flag = 0
			TestZeroh = np.array(np.matrix(eval(str(A).replace("'","").replace(I,str(Kref))))).ravel()
			TestZeroh = sum([0 if i <= 0 else 1 for i in TestZeroh.tolist()])/len(TestZeroh) # Check if some value is <=0
			if TestZeroh < 1: TestZeroh = 0
			TestZerol = np.array(np.matrix(eval(str(A).replace("'","").replace(I,str(1))))).ravel()
			TestZerol = sum([0 if i <= 0 else 1 for i in TestZerol.tolist()])/len(TestZerol) # Check if some value is <=0
			if TestZerol < 1: TestZerol = 0
			while flag == 0 and (TestZeroh+TestZerol) >= 1:
				K = int(Kref)
				print (K)
				TestZero = np.array(np.matrix(eval(str(A).replace("'","").replace(I,str(K))))).ravel()
				TestZero = sum([0 if i <= 0 else 1 for i in TestZero.tolist()])/len(TestZero) # Check if some value is <=0
				if TestZero < 1: TestZero = 0
				GCD = maximumGCD(A, I, K)
				GCD = GCD*TestZero
				if GCD >= K:
					Kref = K
					flag = 1
				else:
					Kref = K-V
			if not Kref:	
				SubsStch = ''
				ProdStch = ''
			else:
				Solution = np.array(np.matrix(eval(str(A).replace("'","").replace(I,str(Kref))))).ravel()
				TestSol = sum([1 for x in Solution if x <=0]) # Check for possible negative values
				if TestSol == 0:
					Solution = [numpy.round(x,1) for x in Solution/min(abs(Solution))]
					LS = len(eq.split('->')[0].split('+'))
					SubsStch = Solution[0:LS]
					ProdStch = Solution[LS:len(Solution)]
				else:
					SubsStch = ''
					ProdStch = ''
		else:
			SubsStch = ''
			ProdStch = ''	         
	else:
		Solution = [numpy.round(x,1) for x in abs(Solution)/min(abs(Solution))]
		LS = len(eq.split('->')[0].split('+'))
		SubsStch = Solution[0:LS]
		ProdStch = Solution[LS:len(Solution)]      
	return SubsStch , ProdStch,eq

## 3.4 New MB_LP   
def MB_LP(eq,AddH,H2O,RxnID):  
	A =  eq2mat(eq).tolist() # Matrix
	C_subs = numpy.ones(shape=(1,len(eq.split("->")[0].split("+")))).tolist()[0]
	C_prod = (numpy.ones(shape=(1,len(eq.split("->")[1].split("+"))))*(-1)).tolist()[0]
	C =C_subs+C_prod # array to minimize
	B = 0
	B = (B,)*len(A) # equaly array
	y=0
	for x in range(0,len(C_subs)):   
		locals()['bound%s' % x] = (1,100) # previous: (-100,0) corrected as the original	          
		y += 1	
	for x in range(y,y+len(C_prod)): 
		locals()['bound%s' % x] = (-100,-1) # previous: (-100,0) corrected as the original      
		y += 1	          
	d=locals()        
	boundaries = [d['bound'+str(x)] for x in range(0,y)]     
	# LP solver
	res = linprog(C, A_eq=A, b_eq=B, bounds=boundaries,options={"disp": False})
	if not res.success or str(res) == "nan":
		SubsStch = ''
		ProdStch = ''
	else:
		Solution = np.round(abs(res.x)/min(abs(res.x)), decimals = 2)
		LS = len(eq.split('->')[0].split('+'))
		SubsStch = Solution[0:LS].tolist()
		ProdStch = Solution[LS:len(Solution)].tolist()       
	return SubsStch , ProdStch,eq

## General MB function:
## Steps:
## 0. Transform Reaction parameters into a string defining the formula of the reaction
## 1. Add missing atoms (Ca,Na,Fe,R,X,K,F)
## 2. Check that all the compounds in the reaction are already defined in the model, otherwise, add
## 3. (Option 1- RxnBalance1) MB reactions using 4 different approache. 
#     The function goes through each method one by one tryin to add H, H2O or both either as substrates or products and stops when a solution is found
## 3. (Option 2- RxnBalance2) MB reactions using 3 different approache. 
##    The function goes through each method one by one tryin to add H, H2O or both either as substrates or products and stops when a solution is found. 
##    If no solution LP approach is used, here all the possible combinations of adding H and/or H2O are calculated and the one with the lowest Stoichimetric indexes is selected
## Option 1: 
def RxnBalance1(eq, RxnID):
	eq_init = eq
	MB = ('', '')
	AddH, AddH2O= '',''
	if len(re.split('\+|->',eq)) < 20:
		eq = AddMissingAtom(eq) # Add atoms in case it is required (Ca,Na,Fe,R,X,K)
		i=0     
#		ListOfFunc = [CountAtom, MB_Core, MB_REM, MB_LP]
		ListOfFunc = [MB_Core, MB_REM, CountAtom, MB_LP]
		AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
		AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it				 
		MB = ListOfFunc[2](eq,AddH,AddH2O,RxnID) #Check if the reaction is already Mass Balanced	
		while not MB[0] and i < len(ListOfFunc): 
			print(ListOfFunc[i])
			AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
			AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it				 
			MB = ListOfFunc[i](eq,AddH,AddH2O,RxnID)			
			print(1)
			if not MB[0] and not re.findall("H[^A-Z0-9]|[^A-Z0-9]H$|^[H][^A-Z0-9]",eq): # not H in the original equatio           
				# Add H to substrates                                   
				AddH = -1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
				AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
				MB = ListOfFunc[i](eq.replace("->","+ H ->"),AddH,AddH2O,RxnID)
				if MB[0]: eq=eq.replace("->","+ H ->")
				print(2)
				# Add H to products
				if not MB[0]:                  
					AddH = 1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
					AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
					MB = ListOfFunc[i](eq.replace("->","-> H +"),AddH,AddH2O,RxnID)
					if MB[0]: eq=eq.replace("->","-> H +")
					print(3)            
			if not MB[0] and not re.findall("H2O[^A-Z0-9]|[^A-Z0-9]H2O$|^[H2O][^A-Z0-9]",eq): # not H2O in the original equation
				# Add H2O to substrates                  
				AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
				AddH2O = -1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
				# Add H2O to products
				MB = ListOfFunc[i](eq.replace("->","+ H2O ->"),AddH,AddH2O,RxnID) 
				if MB[0]: eq=eq.replace("->","+ H2O ->")  
				print(4)         
				if not MB[0]:                  
					AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
					AddH2O = 1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
					MB = ListOfFunc[i](eq.replace("->","-> H2O +"),AddH,AddH2O,RxnID)  
					if MB[0]: eq=eq.replace("->","-> H2O +")
					print(5)
			if not MB[0] and not re.findall("H[^A-Z0-9]|[^A-Z0-9]H$|^[H][^A-Z0-9]",eq) and not re.findall("H2O[^A-Z0-9]|[^A-Z0-9]H2O$|^[H2O][^A-Z0-9]",eq):               
				# Add H2O and H to substrates              
				AddH = -1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
				AddH2O = -1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
				MB = ListOfFunc[i](eq.replace("->","+ H2O + H ->"),AddH,AddH2O,RxnID)
				if MB[0]: eq=eq.replace("->","+ H2O + H ->")  
				print(6)         
				if not MB[0]:                    
					# Add H2O and H to products
					AddH = 1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
					AddH2O = 1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
					MB = ListOfFunc[i](eq.replace("->","-> H2O + H +"),AddH,AddH2O,RxnID)
					if MB[0]: eq=eq.replace("->","-> H2O + H +")
					print(7)
					if not MB[0]:                      
						# Add H2O to products and H to substrates
						AddH = -1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
						AddH2O = 1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
						MB = ListOfFunc[i](eq.replace("->","+ H -> H2O +"),AddH,AddH2O,RxnID) 
						if MB[0]: eq=eq.replace("->","+ H -> H2O +")
						print(8)    
						if not MB[0]:                           
							# Add H2O to substrates and H to products
							AddH = 1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
							AddH2O = -1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
							MB = ListOfFunc[i](eq.replace("->","+ H2O -> H +"),AddH,AddH2O,RxnID)
							if MB[0]: eq=eq.replace("->","+ H2O -> H +")
							print(9)    
			i += 1
#		AllSpecies = ([x for x in [x.replace(" ","") for x in eq.split('->')][0].split('+')],[x for x in [x.replace(" ","") for x in eq.split('->')][1].split('+')])	
		AllSpecies = (([re.sub('^[0-9]*\.?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?','',x.replace(" ","")) for x in [x for x in eq.split('->')][0].split('+')],[re.sub('^[0-9]*','',x.replace(" ","")) for x in [x for x in eq.split('->')][1].split('+')])) 		
		NewSpecies = RxnCompare(eq,eq_init)	      
	if not MB[0]:
		MB = ([1 for x in [x for x in eq.split('->')][0].split('+')],[1 for x in [x for x in eq.split('->')][1].split('+')])
		AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
		AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it                       
		AllSpecies = ([x for x in [x.replace(" ","") for x in eq.split('->')][0].split('+')],[x for x in [x.replace(" ","") for x in eq.split('->')][1].split('+')])
		NewSpecies = ['','']
		eq = eq
	return MB[0], MB[1], AddH, AddH2O, NewSpecies[0], NewSpecies[1], AllSpecies[0], AllSpecies[1], eq



def MissingAtom(eq): 
	Ls = list('abcdefghijklmnopqrstuvwxyz')
	#Mass Balance
	Ss,Es,i,ii=defaultdict(list),[],1,defaultdict(list)
	for p in eq.split('->'):
		for k in p.split('+'):           
			c = [Ls.pop(0), 1]         
			for e,m in re.findall('([A-Z][a-z]?)([0-9]*)',k):              
				m = 1 if m == '' else int(m)
				d = [c[0],c[1]*m*i]
				Ss[e][:0],Es[:0] = [d],[[e,d]]            
				ii[e].append(d[1]) 
		i=-1     
	k = [x[0] for x in list(ii.items()) if x[0] in ['Fe', 'X', 'R', 'Na', 'K', 'Ca', 'F'] and len(x[1]) == 1]
	if k and len(re.split('\+|->',eq)) > 2 and len(re.split('->',eq)[1]) > 2:  
		for x in k:
			if ii[x][0] < 0:
				eq = re.sub(" ->", " + "+x+" ->", eq)
			else:
				eq = re.sub(" ->", " -> "+x+" +", eq)                  
	s= [[x[0],x[1],sum(x[1])] for x in ii.items() if sum(x[1]) != 0]        
	return s



## Option 2:
def mass_balance(eq, RxnID):
	eq_init= eq
	MB = ('', '')
	AddH, AddH2O= '',''
	if len(re.split('\+|->',eq)) < 20:
		eq = AddMissingAtom(eq) # Add atoms in case it is required (Ca,Na,Fe,R,X,K)
		i=0     
#		ListOfFunc = [CountAtom, MB_Core, MB_REM, MB_LP]
		ListOfFunc = [MB_Core, MB_REM, CountAtom, MB_LP]
		AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
		AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it				 
		MB = ListOfFunc[2](eq,AddH,AddH2O,RxnID) #Check if the reaction is already Mass Balanced		
		if MB[0]:
			TestOfBalance = 1
		else:
			TestOfBalance = 0
		while not MB[0] and i < len(ListOfFunc)-1: 
			print(ListOfFunc[i])
			AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
			AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it				 
			MB = ListOfFunc[i](eq,AddH,AddH2O,RxnID)			
			print(1)
			if not MB[0] and not re.findall("H[^A-Z0-9]|[^A-Z0-9]H$|^[H][^A-Z0-9]",eq): # not H in the original equatio           
				# Add H to substrates                                   
				AddH = -1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
				AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
				MB = ListOfFunc[i](eq.replace("->","+ H ->"),AddH,AddH2O,RxnID)
				if MB[0]: eq=eq.replace("->","+ H ->")
				print(2)
				# Add H to products
				if not MB[0]:                  
					AddH = 1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
					AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
					MB = ListOfFunc[i](eq.replace("->","-> H +"),AddH,AddH2O,RxnID)
					if MB[0]: eq=eq.replace("->","-> H +")
					print(3)            
			if not MB[0] and not re.findall("H2O[^A-Z0-9]|[^A-Z0-9]H2O$|^[H2O][^A-Z0-9]",eq): # not H2O in the original equation
				# Add H2O to substrates                  
				AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
				AddH2O = -1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
				# Add H2O to products
				MB = ListOfFunc[i](eq.replace("->","+ H2O ->"),AddH,AddH2O,RxnID) 
				if MB[0]: eq=eq.replace("->","+ H2O ->")  
				print(4)         
				if not MB[0]:                  
					AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
					AddH2O = 1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
					MB = ListOfFunc[i](eq.replace("->","-> H2O +"),AddH,AddH2O,RxnID)  
					if MB[0]: eq=eq.replace("->","-> H2O +")
					print(5)
			if not MB[0] and not re.findall("H[^A-Z0-9]|[^A-Z0-9]H$|^[H][^A-Z0-9]",eq) and not re.findall("H2O[^A-Z0-9]|[^A-Z0-9]H2O$|^[H2O][^A-Z0-9]",eq):               
				# Add H2O and H to substrates              
				AddH = -1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
				AddH2O = -1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
				MB = ListOfFunc[i](eq.replace("->","+ H2O + H ->"),AddH,AddH2O,RxnID)
				if MB[0]: eq=eq.replace("->","+ H2O + H ->")  
				print(6)         
				if not MB[0]:                    
					# Add H2O and H to products
					AddH = 1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
					AddH2O = 1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
					MB = ListOfFunc[i](eq.replace("->","-> H2O + H +"),AddH,AddH2O,RxnID)
					if MB[0]: eq=eq.replace("->","-> H2O + H +")
					print(7)
					if not MB[0]:                      
						# Add H2O to products and H to substrates
						AddH = -1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
						AddH2O = 1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it	
						MB = ListOfFunc[i](eq.replace("->","+ H -> H2O +"),AddH,AddH2O,RxnID) 
						if MB[0]: eq=eq.replace("->","+ H -> H2O +")
						print(8)    
						if not MB[0]:                           
							# Add H2O to substrates and H to products
							AddH = 1 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
							AddH2O = -1 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it
							MB = ListOfFunc[i](eq.replace("->","+ H2O -> H +"),AddH,AddH2O,RxnID)
							if MB[0]: eq=eq.replace("->","+ H2O -> H +")
							print(9)    
			i += 1
#		AllSpecies = ([x for x in [x.replace(" ","") for x in eq.split('->')][0].split('+')],[x for x in [x.replace(" ","") for x in eq.split('->')][1].split('+')])	
#		AllSpecies = (([re.sub('^[0-9]*','',x.replace(" ","")) for x in [x for x in eq.split('->')][0].split('+')],[re.sub('^[0-9]*','',x.replace(" ","")) for x in [x for x in eq.split('->')][1].split('+')])) 		
#		NewSpecies = RxnCompare(eq,eq_init)
		if not MB[0]:
			print(10)    
			i=len(ListOfFunc)-1
			MEtInEq = [x for x in [x.replace(" ","") for x in eq.split('->')][0].split('+')],[x.replace(" ","") for x in [x for x in eq.split('->')][1].split('+')]
			IsH=([1 for x in MEtInEq[0] if x in ['H']],[1 for x in MEtInEq[1] if x in ['H']])
			IsH2O=([1 for x in MEtInEq[0] if x in ['H2O']],[1 for x in MEtInEq[1] if x in ['H2O']])
			MB1 = [ListOfFunc[i](eq,0,0,RxnID),0,0]
			MB1 = [sum(MB1[0][0]+MB1[0][1]),MB1]
			MB2 = MB3 = MB4 = MB5 = MB6 = MB7 = MB8 = MB9 = MB1
			if not IsH[0]:
				MB2 = [ListOfFunc[i](eq.replace("->","+ H ->"),-1,0,RxnID),-1,0]
				MB2 = [sum(MB2[0][0]+MB2[0][1]),MB2]	
			if not IsH[1]:
				MB3 = [ListOfFunc[i](eq.replace("->","-> H +"),1,0,RxnID),1,0]
				MB3 = [sum(MB3[0][0]+MB3[0][1]),MB3]
			if not IsH2O[0]:	
				MB4 = [ListOfFunc[i](eq.replace("->","+ H2O ->"),0,-1,RxnID),0,-1] 
				MB4 = [sum(MB4[0][0]+MB4[0][1]),MB4]
			if not IsH2O[1]:
				MB5 = [ListOfFunc[i](eq.replace("->","-> H2O +"),0,1,RxnID),0,1]  
				MB5 = [sum(MB5[0][0]+MB5[0][1]),MB5]
			if not IsH[0] and not IsH2O[0]:
				MB6 = [ListOfFunc[i](eq.replace("->","+ H2O + H ->"),-1,-1,RxnID),-1,-1] 
				MB6 = [sum(MB6[0][0]+MB6[0][1]),MB6]
			if not IsH[1] and not IsH2O[1]:
				MB7 = [ListOfFunc[i](eq.replace("->","-> H2O + H +"),1,1,RxnID),1,1]
				MB7 = [sum(MB7[0][0]+MB7[0][1]),MB7]
			if not IsH[0] and not IsH2O[1]:	
				MB8 = [ListOfFunc[i](eq.replace("->","+ H -> H2O +"),-1,1,RxnID),-1,1]  
				MB8 = [sum(MB8[0][0]+MB8[0][1]),MB8]
			if not IsH[1] and not IsH2O[0]:			
				MB9 = [ListOfFunc[i](eq.replace("->","+ H2O -> H +"),1,-1,RxnID),1,-1]
				MB9 = [sum(MB9[0][0]+MB9[0][1]),MB9]	
			MBOpt = sorted([MB1,MB2,MB3,MB4,MB5,MB6,MB7,MB8,MB9])[0][1]
			MB = MBOpt[0]
			AddH = MBOpt[1]
			AddH2O = MBOpt[2]
			eq = MB[2]
#		AllSpecies = ([x for x in [x.replace(" ","") for x in eq.split('->')][0].split('+')],[x for x in [x.replace(" ","") for x in eq.split('->')][1].split('+')])
		AllSpecies = (([re.sub('^[0-9]*\.?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?','',x.replace(" ","")) for x in [x for x in eq.split('->')][0].split('+')],[re.sub('^[0-9]*','',x.replace(" ","")) for x in [x for x in eq.split('->')][1].split('+')])) 		
		NewSpecies = RxnCompare(eq,eq_init)
	s=''        
	if not MB[0]:
		MB = ([1 for x in [x for x in eq.split('->')][0].split('+')],[1 for x in [x for x in eq.split('->')][1].split('+')])
		AddH = 0 # In case it is necessary to add H+, and H+ compound has not been previously added into the model, this flag is used to include it
		AddH2O = 0 # In case it is necessary to add H2O, and H2O compound has not been previously added into the model, this flag is used to include it                       
#		AllSpecies = ([x for x in [x.replace(" ","") for x in eq.split('->')][0].split('+')],[x for x in [x.replace(" ","") for x in eq.split('->')][1].split('+')])
		AllSpecies = (([re.sub('^[0-9]*\.?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?[0-9]?','',x.replace(" ","")) for x in [x for x in eq.split('->')][0].split('+')],[re.sub('^[0-9]*','',x.replace(" ","")) for x in [x for x in eq.split('->')][1].split('+')])) 		
		NewSpecies = ['','']
		TestOfBalance = 2
	return MB[0], MB[1], AddH, AddH2O, NewSpecies[0], NewSpecies[1], AllSpecies[0], AllSpecies[1], eq, eq_init, TestOfBalance