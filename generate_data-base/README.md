# generate_database

2. generate_database: generates a metabolic network based using current available information from a large number of molecular databases:

	- generate_database.py:
		- Mass balance metabolic reactions (equations_generate_database.mass_balance)
		- Constritruction of SGPRs and GPRs (gpr.auth_gpr.getGPR)
		- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)