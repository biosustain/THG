# Build model

## Brief Protocol Description 

3. build model: Used to enrich gene annotation, curate and expands a metabolic network model

	- build_model.py:
		- Mass balance metabolic reactions (equations_build_model.mass_balance)
		- Constritruction of SGPRs and GPRs (gpr.auth_gpr.getGPR)
		- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
		- Other functions for build_model.py (function_build_model.py)

## Usage

To run the model building, use the following command:
		
The code can be run as:

```
python3 build_model.py
```
