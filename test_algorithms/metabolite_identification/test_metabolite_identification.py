"""Test metabolite annotation."""

import pytest
import pandas as pd
from os.path import join, dirname
from functions_metabolite_identification import identify_metabolite


PUBCHEM_FIXTURE = list(
    pd.read_csv(join("files", "kegg.tsv"), sep="\t").itertuples(
        index=False, name=None
    )
)

@pytest.mark.parametrize("name,formula,identifier,expected", PUBCHEM_FIXTURE)
def test_all_metabolites_can_be_identified(name:str, formula: str, identifier: str, expected: str):
    """Check that there was a match against pubchem and identifier is right."""
    met_result = identify_metabolite(name, formula, identifier, 0.70)
    assert met_result is not None
    assert str(expected) in met_result
