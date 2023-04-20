# Build model

This repository is a python module to enrich gene annotation, curate and expands a metabolic network model based on information currently available in multiple online databases.

## Usage

To run the script, use the following command:

```
python3 build_model.py
```


3. Build model: 

	- build_model.py:
		- Mass balance metabolic reactions (equations_build_model.mass_balance)
		- Constritruction of SGPRs and GPRs (gpr.auth_gpr.getGPR)
		- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
		- Other functions for build_model.py (function_build_model.py)

```
