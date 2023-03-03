import difflib
import logging
import re

import cobra
import pandas as pd
import pubchempy as pcp
import numpy as np
import pickle
from collections import defaultdict
from typing import List, Dict, Optional, Tuple
from function_reac_identification import *


LOGGER = logging.getLogger(__name__)
H_PATTERN = re.compile(r"H[0-9]+")
CHARGE_PATTERN = re.compile(r"[-\+][0-9]+")

def global_met_annotation_file():
    global annotation_file
    annotation_file = 'reports/met_annotation.tsv'
    return annotation_file


def gather_metabolites(
    model: cobra.Model) -> List[Tuple[str, str, str, str]]:
    """Gather metabolites info from a `model`, filtering out from `identifiers`."""
    met_list = []
    for met in model.metabolites:
        met_tuple = (met.name, met.formula, met.annotation, met.id[:-1])
        if met_tuple not in met_list:
            met_list.append(met_tuple)
    return met_list


def formula_similarity(formula_left: str, formula_right: str) -> float:
    """Compute similiratiy based on the stdlib SequenceMatcher."""
    pruned_left = H_PATTERN.sub("H", CHARGE_PATTERN.sub("", atom(formula_left)))
    pruned_right = H_PATTERN.sub("H", CHARGE_PATTERN.sub("", atom(formula_right)))
    return difflib.SequenceMatcher(
        lambda x: x == " ", pruned_left, pruned_right
    ).ratio()


def identify_metabolite(
    name: str, formula: str, iden: str, threshold: float = 0.82
) -> Optional[str]:
    """Try to match metabolite info to a pubchem compound to get the annotation.

    Parameters
    ----------
    name: str
    formula: str
    iden: str
    threshold: float, default=0.82
        minimum similarity ratio between the molecular formulas. The default (0.82)
        was decided by performing sensibility tests of the metabolites in Human 1.
    
    Returns
    -------
    met_result: Optional[str]
        Tab-separated str. If None, the metabolite was not identified.

    """
    CompoundID = pcp.get_cids(name.strip(), "name")
    met_result = None
    if not CompoundID:
        LOGGER.warning("get_cids did not work")
        CompoundID = pcp.get_cids(re.sub(r"[\(\)]", "", name), "name")
    if CompoundID:
        Compound = pcp.Compound.from_cid(CompoundID)
        molecular_formula = Compound.molecular_formula

        if formula_similarity(formula, molecular_formula) >= threshold:
            Compound = pcp.Compound.from_cid(CompoundID)
            CompoundID = Compound.cid
            synonyms = Compound.synonyms
            Metabolite = sorted(
                [
                    [
                        difflib.SequenceMatcher(
                            None, name.lower(), synonym.lower()
                        ).ratio(),
                        synonym,
                    ]
                    for synonym in synonyms
                ]
            )[-1][1]
            Kegg = "".join(list(filter(re.compile("[CG][0-9]{5}$").match, synonyms)))
            LIPIDMAPSID = "".join(
                list(filter(re.compile("L[A-Z]{3}[0-9]+$").match, synonyms))
            )
            CHEBI = "".join(list(filter(re.compile("CHEBI:[0-9]+$").match, synonyms)))
            inchi = Compound.inchi
            inchikey = Compound.inchikey
            met_result = (
                Metabolite
                + "\t"
                + ""
                + "\t"
                + molecular_formula
                + "\t"
                + LIPIDMAPSID
                + "\t"
                + Kegg
                + "\t"
                + CHEBI
                + "\t"
                + str(CompoundID)
                + "\t"
                + ""
                + "\t"
                + inchikey
                + "\t"
                + inchi
                + "\t"
                + ""
                + "\t"
                + iden
                + "\n"
            )
    else:
        LOGGER.warning("get_cids extra did not work")
    return met_result

    
def generate_met_annotation(
    met_list: List[Tuple[str, str, str, str]], out: str = global_met_annotation_file()
) -> Tuple[List, List]:
    """Identify and write metabolite annotations to file (Tab-separated).

    Returns
    -------
    unnanotated: list[tuple[str, str, str, str]]
        list of unnanotated metabolites
    anotated: list[tuple[str, str, str, str]]
        list of anotated metabolites
    """
    unnanotated, annotated = [], []

    with open(out, "w") as f:
        for name, formula, annotation, iden in met_list:
            try:
                result = identify_metabolite(name, formula, iden)
                met = [name, formula, annotation, iden]
                if result is None:
                    unnanotated.append(met)
                else:
                    f.write(result)
                    annotated.append(met)
            except Exception as e:
                print(e)
                continue
    return annotated, unnanotated


def remove_null_value(d):
    return {
        a: c
        for a, b in d.items()
        if (c := (b if not isinstance(b, dict) else remove_null_value(b)))
    }


def process_annotation(
    annotation_file: str = global_met_annotation_file()) -> Dict:
    """Generate metabolite annotation file.

    Parameters
    ----------
    annotation_file: str
        Tab-separated file, generated with `generate_met_annotation`
    """
    variableFile = pd.read_csv(annotation_file, sep="\t", header=None)
    variableFile[5] = variableFile[5].str.extract("(CHEBI:[0-9]+)", expand=True)
    variableFile = variableFile.rename(
        columns={
            3: "lipidmaps",
            4: "kegg.compound",
            5: "chebi",
            6: "pubchem.compound",
            8: "inchikey",
            9: "inchi",
        }
    )
    variableFile = variableFile.fillna(str())

    annotation = defaultdict(dict)

    for _, row in variableFile.iterrows():  
        annotation[row[11]] = row.drop([0, 1, 2, 7, 10, 11]).to_dict()
    annotation = remove_null_value(annotation)
    return annotation


def atom(Formula):
    try:
        if "C" in Formula:
            C = re.findall(r"(C[0-9]+)", Formula)
            if not C:
                C = ["C1"]
        else:
            C = []
        if "H" in Formula:
            H = re.findall(r"(H\+*[0-9]+)", Formula)
            if not H:
                H = ["H1"]
        else:
            H = []
        if "O" in Formula:
            O = re.findall(r"(O[0-9]+)", Formula)
            if not O:
                O = ["O1"]
        else:
            O = []
        if "N" in Formula:
            N = re.findall(r"(N[0-9]+)", Formula)
            if not N:
                N = ["N1"]
        else:
            N = []
        if "P" in Formula:
            P = re.findall(r"(P[0-9]+)", Formula)
            if not P:
                P = ["P1"]
        else:
            P = []
        if "S" in Formula:
            S = re.findall(r"(Se?[0-9]+)", Formula)
            if not S:
                S = ["S1"]
        else:
            S = []
        if "I" in Formula:
            I = re.findall(r"(I[0-9]+)", Formula)
            if not I:
                I = ["I1"]
        else:
            I = []
        if "F" in Formula:
            F = re.findall(r"(Fe*[0-9]+)", Formula)
            if not F:
                F = ["F1"]
        else:
            F = []
        if "R" in Formula:
            R = re.findall(r"(R[0-9]+)", Formula)
            if not R:
                R = ["R1"]
        else:
            R = []
        return "".join(str(i) for i in C + H + O + N + P + S + I + F + R)
    except Exception as e:
        print(e)