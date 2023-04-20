# Build model

The build model repository is a python module to enrich gene annotation, curate and expands a metabolic network model based on information currently available in multiple online databases.

```build_model.py``` is the main script and will execute the following tasks when run:
- Mass balance metabolic reactions (equations_build_model.mass_balance)
- Construct SGPRs and GPRs (gpr.auth_gpr.getGPR)
- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
- Other tasks (function_build_model.py)

## Usage

To run the script, use the following command:

```
python3 build_model.py
```
