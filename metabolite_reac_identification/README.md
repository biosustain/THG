# metabolite_reac_identification

1. metabolite_reac_identification: Enrich metabolite and reaction annotation of a reference model
	- metabolite_reac_identification.py:
		- Identifies metabolites by name and molecular formula (function_metabolite_identification.generate_met_annotation)
		- Identifies reactions by unique combination of substrates and products (function_reac_identification.execute_jaccard)
		- Writes the new model (function_annotate_cobra_model.annotate_cobra_model)
