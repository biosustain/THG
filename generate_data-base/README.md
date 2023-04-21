# Generate databae

The generate database repository is a python module to generates a metabolic network based using current available information from a large number of online databases.

```generate_database.py``` is the main script and will execute the following tasks when run:
- Mass balance metabolic reactions (equations_generate_database.mass_balance)
- Construct SGPRs and GPRs (gpr.auth_gpr.getGPR)
- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
- Other tasks (functions_generate_database.py)
- Other tasks (class_generate_database.py)

## Usage

To run the script, use the following command:

```
python3 generate_database.py
```

## Models and files

The models and files can be found in the respective folder where a brief description can be read.  
