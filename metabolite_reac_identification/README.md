# metabolite_reac_identification

1. metabolite_reac_identification: Enrich metabolite and reaction annotation of a reference model
	- metabolite_reac_identification.py:
		- Identifies metabolites by name and molecular formula (generate_met_annotation)
		- Identifies reactions by unique combination of substrates and products (execute_jaccard)
	- function_metabolite_identification.py: 
		- The functions to identify metabolites by name and molecular formula (generate_met_annotation)
	- function_reac_identification.py:
		- The functions to identifify reactions by unique combination of substrates and products (execute_jaccard)
	- function_annotate_cobra_model.py:
		- The function to write the new model

2. models:
	- Human-GEM_2022-06-21.xml: Used reference model
	- Human Database.xml:Used reference database (developed in THG/generate_data-base)
	- THG-beta1.1.xml: New model with enriched annotation (version 1.1)	
	- THG-beta1.1.1.xml: New model with enriched annotation (version 1.1.1)
