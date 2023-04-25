# Metabolite reac identification

The "metabolite_reac_identification" repository is a python module to enrich metabolite and reaction annotation of a reference human genome-scale metabolic model

metabolite_reac_identification.py:
- Identifies metabolites by name and molecular formula (function_metabolite_identification.generate_met_annotation)
- Identifies reactions by unique combination of substrates and products (function_reac_identification.execute_jaccard)
- Correct/Enrich reference model with the information collected in the previous steps
- Writes the corrected/enriched new model (function_annotate_cobra_model.annotate_cobra_model)


## Running metabolite reac identification

To run the metabolite reac identification:

```
python3 metabolite_reac_identification.py
```

