# merge metabolic network and network consistency

4. merge metabolic network and network consistency: merge two metabolic network while ensuring the consistency of the new model in a interative fashion

	- merge_metabolic_networks.py:
		- Mass balance metabolic reactions (functions_mass_balance.mass_balance)
		- MEMOTE function that checks is a GEM is stoichiometricaly consistent (functions_network_consistency.test_stoichiometric_consistency)
		- Other functions for merge_metabolic_networks.py (functions_merge_metabolic_networks.py) 
