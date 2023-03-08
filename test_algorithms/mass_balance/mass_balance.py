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
from cobra.io import read_sbml_model,write_sbml_model
#from math import gcd
import math
import csv
import pandas as pd
from equations_mass_balance import *

cobra_model = read_sbml_model('models/F-H1-compartments_DB-2022-10-05.xml')
writer = pd.ExcelWriter('files/mass_balance.xlsx', engine='xlsxwriter')

ListOfMetFrom = {k: str() if v is None else v for (k, v) in {x.id: x.formula for x in cobra_model.metabolites}.items()}
v = {8: 'H1', 11: 'H1', 13: 'DB'}

RxnList, Rxn2MB = [],[]
for x in cobra_model.reactions:    
	eq = re.sub('(^| )[0-9\.]+','',re.sub('e-[0-9]+','', re.sub('<=>','->', re.sub('-->','->', re.sub('  ',' ', re.sub('[0-9\.]+e-[0-9]+','', x.reaction)))))).strip() # correct 
	Rxn = re.sub('[a-z]+[0-9]*','', eq) # correct
	eq = re.sub('<=>', '->', re.sub('-->','->', re.sub('e-[0-9]+','', re.sub('(^| )[0-9\.]+','', x.reaction).strip()).strip() ) )   
	Rxn = re.sub('[a-z]+[0-9]*','', re.sub('e-[0-9]+','', re.sub('(^| )[0-9\.]+','', x.reaction).strip()).strip() )         
	RxnID =  x.annotation.get('kegg.reaction',str())        
	species = [x.id for x in x.reactants] + [x.id for x in x.products]     
	if not Rxn in RxnList:   
		for y in species:             
			eq = re.sub(y, ListOfMetFrom.get(y,str()), eq)                 
		Rxn2MB.append([x.id, RxnID, v[len(x.id)], Rxn, eq])     
		RxnList.append(Rxn) 
        
## Apply the MB algorithm

AddH = 0
H2O = 0
RxnID = 'R'
MBSummary = []
MB_e = []

for i in range(0, len(Rxn2MB)):
	eq = Rxn2MB[i][-1]
	y = Rxn2MB[i][0]
	try:                
		b = mass_balance(eq,RxnID)   
		MBSummary.append([b[0],b[1],b[6],b[7],b[9],b[10]])      
	except Exception as e:          
		reactants = [ListOfMetFrom.get(x.id, str()) for x in cobra_model.reactions.get_by_id(y).reactants]
		products = [ListOfMetFrom.get(x.id, str()) for x in cobra_model.reactions.get_by_id(y).products]        
		MB_e.append([RxnID, eq])   
		MBSummary.append([[],[],reactants,products,eq,3])   
        
        
## Write the results to file 

MBSum = pd.DataFrame(MBSummary)
MBSum.columns = [5,6,7,8,4,9]
#MBSum = pd.merge(pd.DataFrame(Rxn2MB), MBSum, on=[4]).drop_duplicates([3]) # correct 
MBSum = pd.merge(pd.DataFrame(Rxn2MB), MBSum, on=[4]).drop_duplicates([4])
MBSum = MBSum.astype(str)
MBSum[(MBSum[9].str.contains('0',na=False)) | (MBSum[9].str.contains('1',na=False))].drop([9],axis=1).to_excel(writer, sheet_name='mass_balanced', index=False)
MBSum[(MBSum[9].str.contains('2',na=False)) | (MBSum[9].str.contains('3',na=False))].drop([9],axis=1).to_excel(writer, sheet_name='not_mass_balanced', index=False)
writer.close()