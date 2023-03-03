import re
import pandas as pd
from io import StringIO
import pickle
from typing import Optional, Union


def process_reac(input_file, MetID, MetIDH, MetIDH2O: Optional = '', take=-1):
    result = ""
    variableFile = open(input_file, "r").read()
    for n, reaction in enumerate(re.findall(r"<reaction.+?<\/reaction>", variableFile, re.DOTALL)):
        if take < n:
            result += identify_reaction(reaction, n, MetID, MetIDH, MetIDH2O)
    RxnJI = StringIO(result)
    out = pd.read_csv(RxnJI, sep=" ", header=None)
    out.dropna(subset=[5,6],inplace=True)
    return out


def gather_kegg_metabolites(Output):
    Metabolite = {}
    with open(Output) as m:
        m = m.read()
        for species in re.findall("<species metaid.+?<\/species>", m, re.DOTALL):
            CompoundID = ""
            MetID = re.findall('about="#(.+?)">', species)[0]
            if re.findall("compound.([A-Z][0-9]+)", species):
                CompoundID = re.findall("compound.([A-Z][0-9]+)", species)[0]
            if MetID and CompoundID and not MetID in Metabolite:
                Metabolite[MetID] = CompoundID + re.findall("[a-z]+[0-9]*", MetID)[0]
    return Metabolite


def replace_met_id_by_met_kegg(
    infile, dictionary):
    infile = infile[infile[5].notna() & infile[6].notna()]
    infile[5].replace(dictionary, regex=True, inplace=True)
    infile[6].replace(dictionary, regex=True, inplace=True)
    return infile


def execute_jaccard(b,c2): 
    d = ''
    for l,r,a in zip(b[5],b[6],b[0]):
        for l2,r2,b in zip(c2[5],c2[6],c2[1]):
            a, b = str(a), str(b)
            Jllrr = jaccard(l.split(','),l2.split(','))+jaccard(r.split(','),r2.split(','))
            Jlrrl = jaccard(l.split(','),r2.split(','))+jaccard(r.split(','),l2.split(','))
            if Jllrr == 2.0 or Jlrrl == 2.0:
                d = d + f"{a}\t{l}\t{r}\t{b}\t{l2}\t{r2}\t2.0\n"
                break
    d = StringIO(d)
    return pd.read_csv(d, sep='\t', header=None)


def process_jaccard(left_model, right_model, jaccard_result):
    
    """Postprocess dataframe from ReactionJI."""
    """Write reaction annotation file as dict """
    right_model.columns= [13,14,15,16,17,18,19,20,21,22,23,24]
    jaccard_result.columns = [0,5,6,14,18,19,12]
    jaccard_result = pd.merge(pd.merge(jaccard_result,right_model,on=[18,19,14]), left_model, on=[0,5,6]).drop(columns=[11])
    jaccard_result = jaccard_result.sort_index(axis=1)
    jaccard_result = jaccard_result.sort_values(by=[0,1])
    jaccard_result[100] = jaccard_result[1]+jaccard_result[2]
    mm = jaccard_result[(jaccard_result[100].isin(jaccard_result[100][jaccard_result[100].duplicated()]))]
    mm = mm[(mm[8] == mm[21])]
    mm = mm[(mm[100].isin(mm[100][mm[100].duplicated()]))]
    mm = mm.drop_duplicates(100).reset_index(drop=True) 
    jaccard_result = jaccard_result.loc[~jaccard_result[100].isin(mm[100])]
    jaccard_result = jaccard_result[(jaccard_result[8] == jaccard_result[21])]
    jaccard_result=jaccard_result[jaccard_result[1]!=jaccard_result[14]]
    jaccard_result = jaccard_result.rename(columns={1: "kegg.reaction"})
    jaccard_result = jaccard_result[['kegg.reaction',14, 24]]
    r = dict()
    for i, row in jaccard_result.iterrows():
        r[row[24]] = row.drop([14,24]).to_dict() 
    return r


def identify_reaction(reaction: str, n: int, MetID, MetIDH, MetIDH2O: Optional = ''):
    MARID = "MAR[0-9]+"
    MetID2 = "([A-Z0-9]+)"
    s = "<listOfReactants>(.+?)<.listOfReactants>"
    p = "<listOfProducts>(.+?)<.listOfProducts>"
    p2, c2 = [], []
    RxnID, MAR, e = "", "", ""
    if re.findall("/(R[0-9]+)", reaction):
        RxnID = re.findall("/(R[0-9]+)", reaction)[0]
    if re.findall(MARID, reaction):
        MAR = re.findall(MARID, reaction)[0]
    s2 = re.findall(s, reaction, re.DOTALL)[0]
    s2 = re.findall(MetID, s2)
    c1 = list(set([re.findall("[a-z]+[0-9]*", x)[0] for x in s2]))
    s2 = [
        x
        for x in s2
        if re.findall(MetID2, x)[0] not in MetIDH
        and re.findall(MetID2, x)[0] not in MetIDH2O
    ]
    if re.findall(p, reaction, re.DOTALL):
        p2 = re.findall(p, reaction, re.DOTALL)[0]
        p2 = re.findall(MetID, p2)
        c2 = list(set([re.findall("[a-z]+[0-9]*", x)[0] for x in p2]))
        p2 = [
            x
            for x in p2
            if re.findall(MetID2, x)[0] not in MetIDH
            and re.findall(MetID2, x)[0] not in MetIDH2O
        ]
    sp = (",".join(map(str, s2)) + " " + ",".join(map(str, p2))).replace("_", "")
    c = ",".join(list(set(re.findall("[a-z]+[0-9]*", sp))))
    c2 = ",".join(list(set(c1 + c2)))
    Stch = [
        str(round(float(re.findall(x + '" stoichiometry="(.+?)"', reaction)[0])))
        for x in s2 + p2
    ]
    Stch = ",".join(Stch)
    Stch2 = [
        str(round(float(x))) for x in re.findall('stoichiometry="(.+?)"', reaction)
    ]
    Stch2 = ",".join(Stch2)
    Rev = re.findall('reversible="(.+?)"', reaction)[0]
    e = ",".join(re.findall(r'\/[a-z]+-[a-z]+\/\+?([\.0-9]+)"/>', reaction))
    g = ",".join(re.findall('geneProduct="(.+?)"/>', reaction))
    return (
        str(n)
        + " "
        + RxnID
        + " "
        + c
        + " "
        + c2
        + " "
        + Rev
        + " "
        + sp
        + " "
        + Stch
        + " "
        + Stch2
        + " "
        + e
        + " "
        + g
        + " "
        + MAR
        + "\n"
    )


def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union