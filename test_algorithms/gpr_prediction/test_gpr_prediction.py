"""Test that GPR are produced from EC-Numbers categorized by EC-number."""

import re
from os.path import dirname, join

import pandas as pd
import pytest
import requests

from functions_auth_gpr import getGPR

OP_PAT = re.compile("and|or")
GPR_FIXTURES = [
    list(
        pd.read_csv(
            join("files", f"gprs_ec{i}.tsv"), sep="\t"
        ).itertuples(index=False, name=None)
    )
    for i in range(1, 10)
    # EC7 not in data
    if i != 7
]


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[0])
def test_gprs_ec1(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 1."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[1])
def test_gprs_ec2(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 2."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[2])
def test_gprs_ec3(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 3."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[3])
def test_gprs_ec4(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 4."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[4])
def test_gprs_ec5(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 5."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[5])
def test_gprs_ec6(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 6."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[6])
def test_gprs_ec8(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 8."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr


@pytest.mark.parametrize("ecnumber,expected_gprs", GPR_FIXTURES[7])
def test_gprs_ec9(session: requests.Session, ecnumber: str, expected_gprs: str):
    """Check expected GPRs in group of EC class 9."""
    fetched_gpr = getGPR(ecnumber, session)
    assert fetched_gpr is not None
    genes, gene_names, gene_ids, sgpr, gpr = fetched_gpr
    for expected_gpr in OP_PAT.split(expected_gprs):
        assert expected_gpr.strip().replace("(", "").replace(")", "") in gpr
