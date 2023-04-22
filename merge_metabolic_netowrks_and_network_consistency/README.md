# Merge metabolic network and network consistency

The merge metabolic networks and network consistency repository is a python module to merge two metabolic network while ensuring the consistency of the new model in a interative fashion.

merge_metabolic_networks.py:
- Mass balance metabolic reactions (functions_mass_balance.mass_balance)
- MEMOTE function that checks is a GEM is stoichiometricaly consistent (functions_network_consistency.test_stoichiometric_consistency)
- Additional functions necessary to run the script (functions_merge_metabolic_networks.py)


## Running merge metabolic network and network consistency

To run the metabolic network merging and network consistency:

```
python3 merge_metabolic_networks.py
```
