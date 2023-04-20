# Metabolite reac identification

This repository is a python module to enrich metabolite and reaction annotation of a reference model

```metabolite_reac_identification.py``` is the main script and will execute the following tasks when run:
- Identifies metabolites by name and molecular formula (function_metabolite_identification.generate_met_annotation)
- Identifies reactions by unique combination of substrates and products (function_reac_identification.execute_jaccard)
- Writes the new model (function_annotate_cobra_model.annotate_cobra_model)


## Usage

To run the script, use the following command:

```
metabolite_reac_identification.py
```

## Models and files

The models and files can be found in the respective folder where a brief description can be read.
