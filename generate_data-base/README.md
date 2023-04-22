# Generate databae

The generate database directory is a python module to generates a metabolic network based using current available information from a large number of online databases.

```generate_database.py```:
- Mass balance metabolic reactions (equations_generate_database.mass_balance)
- Construct SGPRs and GPRs (gpr.auth_gpr.getGPR)
- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
- Additional functions necessary to run the script (functions_generate_database.py)
- Classes necessary to run the script (class_generate_database.py)

## Usage

To run the script, use the following command:

```
python3 generate_database.py
```
