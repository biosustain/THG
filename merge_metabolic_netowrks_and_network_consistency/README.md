# Merge metabolic network and network consistency

This repository is a python module to merge two metabolic network while ensuring the consistency of the new model in a interative fashion.

```merge_metabolic_networks.py``` is the main script and will use the following modules when run:
- Mass balance metabolic reactions (functions_mass_balance.mass_balance)
- MEMOTE function that checks is a GEM is stoichiometricaly consistent (functions_network_consistency.test_stoichiometric_consistency)
- Other functions for merge_metabolic_networks.py (functions_merge_metabolic_networks.py, functions_network_consistency.py) 


## Usage

To run the script, use the following command:

```
merge_metabolic_networks.py
```
