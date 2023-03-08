"""Test that annotation of reactions works as expected."""

from os.path import dirname, join
from typing import Dict

import cobra
import pandas as pd
import pytest

from functions_reac_identification import (
    execute_jaccard,
    gather_kegg_metabolites,
    process_reac,
    identify_reaction,
    process_jaccard,
)

MetIDH = "MAM02040"
MetIDH2O = "MAM02039"
MetID = "[A-Z]+[0-9]+[a-z]+[0-9]*"

ECNUMBER_FIXTURE = list(
    pd.read_csv(join("files", "ec-number.tsv"), sep="\t").itertuples(
        index=False, name=None
    )
)
# the reaction map was parsed from the 00_Lipid_Carma_KEGG_Model.xml, filtered by
# those who had ECNUMBER and then deduplicated by lipid identifier. Not all reactions
# are in the result of the Jaccard identification
REACMAP_FIXTURE = list(
    pd.read_csv(join("files", "reaction_map.csv")).itertuples(
        index=False, name=None
    )
)


@pytest.mark.parametrize("sbml_reaction,expected", ECNUMBER_FIXTURE)
def test_all_reactions_return_right_ecnumber(sbml_reaction: str, expected: str):
    """Check that the EC number is found."""
    assert expected in identify_reaction(sbml_reaction, 0, MetID, MetIDH2O, MetID)


def test_reaction_jaccard_works(
    human_path: str, lipids_path: str, ground_truth_reacs_small: str
):
    """Check that the expected matches are found, EC number matching."""
    ground: pd.DataFrame = pd.read_csv(ground_truth_reacs_small)
    MetIDH = "MAM02040"
    MetIDH2O = "MAM02039"
    MetID = "[A-Z]+[0-9]+[a-z]+[0-9]*"
    human = ReactionJI(human_path, MetID, MetIDH, MetIDH2O).sort_values(by=[5, 6])
    lipids = ReactionJI(
        lipids_path, r"([A-Z][0-9]+_?[a-z]+[0-9]*)", "C00080", "C00001"
    ).sort_values(by=[5, 6])
    MetID = MetID2(human_path)
    human = human[human[5].notna() & human[6].notna()]
    human[5], human[6] = human[5].replace(MetID, regex=True), human[6].replace(
        MetID, regex=True
    )
    jaccard_result = JI(lipids, human)
    processed_result = process_jaccard(lipids, human, jaccard_result)
    model = cobra.io.read_sbml_model(human_path)
    # the EC-NUMBER comes fromt 00_Lipid
    assert all(
        ground.apply(
            lambda x: x.loc["ec-number"]
            in model.reactions.get_by_id(processed_result[x.lipids][0]).annotation[
                "ec-code"
            ],
            axis=1,
        )
    )


@pytest.mark.parametrize("reac,expected_ec", REACMAP_FIXTURE)
def test_reactions_are_the_same(
    human_model: cobra.Model, processed_jaccard: Dict, reac: str, expected_ec: str
):
    """Check that the expected matches are found, EC number matching.

    The processed_jaccard was precomputed with the `process_jaccard()` function.
    """
    for comp in {'c', 'm', 'n', 'x', 'e', 'g', 'l', 'r'}:
        reac_id = reac + comp
        if reac_id in processed_jaccard:
            break
    else:
        pytest.skip(f"{reac} is not part of the Jaccard result.")
    if "ec-code" not in human_model.reactions.get_by_id(processed_jaccard[reac_id][0]).annotation:
        pytest.skip(f"{reac} cannot be tested (not ec number in human model).")
    assert (
        expected_ec
        in human_model.reactions.get_by_id(processed_jaccard[reac_id][0]).annotation[
            "ec-code"
        ]
    )
