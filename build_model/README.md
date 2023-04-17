# build model

3. build model: Enrich gene annotation, curate and expands a metabolic network model

	- build_model.py:
		- Mass balance metabolic reactions (equations_build_model.mass_balance)
		- Constritruction of SGPRs and GPRs (gpr.auth_gpr.getGPR)
		- Identifies the cellular location of the metabolic reactions (gpr.getLocation.getLocation)
