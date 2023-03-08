"""Top-level fixtures for testing."""

import json
from os.path import dirname, join
from typing import Dict

import cobra
import pytest


@pytest.fixture(scope="module")
def human_path() -> str:
    """Path to human GEM."""
    return join(dirname(__file__), "data", "Human-GEM.xml")


@pytest.fixture
def lipids_path() -> str:
    return join(dirname(__file__), "data", "00_Lipid_Carma_KEGG_Model.xml")


@pytest.fixture
def ground_truth_reacs_small() -> str:
    return join(dirname(__file__), "data", "reaction_map.csv")


@pytest.fixture(scope="module")
def processed_jaccard() -> Dict:
    with open(join(dirname(__file__), "data", "processed.json")) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def human_model(human_path: str) -> cobra.Model:
    return cobra.io.read_sbml_model(human_path)
