# Generate_database

## Brief Protocol Description 		


2. generate_database: The script generates a metabolic network based using current available information from a large number of online databases:

	- generate_database.py:
		- Mass balance metabolic reactions (equations_generate_database.mass_balance)
		- Constritruction of SGPRs and GPRs (gpr.auth_gpr.getGPR)
		- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
		- Other functions for generate_database.py (functions_generate_database.py)
		- Classes for generate_database.py (class_generate_database.py)

## Usage

To run the generate database script, use the following command:

```
python3 generate_database.py
```
