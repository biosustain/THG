# THG
A protocol for the automatic construction of highly curated genome-scale models of human metabolism.

The protocol enables the automatic curation and/or expansion of existing human GEMs or generates a highly curated metabolic network based on the current information retrieved from multiple databases in real time.

The different steps in the protocol can be found in different folders and can be executed independently


1. metabolite_reac_identification: Enrich gene, metabolite and reaction annotation of a reference model
	- metabolite_reac_identification.py: 
		- Identifies metabolites by name and molecular formula (generate_met_annotation)
		- Identifies reactions by unique combination of substrates and products (execute_jaccard)

2. generate_database: generates a metabolic network based using current available information from a large number of molecular databases:
	- generate_database.py:
		- Mass balance metabolic reactions (mass_balance)
		- Constritruction of SGPRs and GPRs (getGPR)
		- Identifies the cellular location of the metabolic reactions (getLocation)

3. build model: curate and expands a metabolic network model
	- build_model.py:
		- Mass balance metabolic reactions (mass_balance)
		- Constritruction of SGPRs and GPRs (getGPR)
		- Identifies the cellular location of the metabolic reactions (getLocation)

4. merge metabolic network and network consistency: merge two metabolic network while ensuring the consistency of the new model in a interative fashion
	- merge_metabolic_networks.py:
		- Mass balance metabolic reactions (mass_balance)
		- Constritruction of SGPRs and GPRs (getGPR)
		- Identifies the cellular location of the metabolic reactions (getLocation)

		- MEMOTE function that checks is a GEM is stoichiometricaly consistent (test_stoichiometric_consistency)

5. memote_and_task_analysis: performs a MEMOTE analysis to assess the quality of the metabolic model and evaluates if the the model can perform certain cellular functions (tasks)
	- see README.md in memote_and_task_analysis folder

6. compare_models: genereates a report comparing two GEMs:
	- compare_models.py

7. generate_figures: combines information from GEMs SBML files and the reports generated in "compare_models" to generate figures
	- create_figures.py








# get_gpr and get_location

First, install the repository:

```bash
pip install .
```

After that, you need to create a BioCyc account. This works if you are in an institution that has BioCyc licenses, usually requiring physical presence.

Once installed and with a BioCyc account, the following script will fetch different ec numbers in a list:

```python
from gpr.auth_gpr import getGPR, setup_biocyc_session

ec_numbers = ["1.1.1.1", "1.1.1.204"]
# this will make authenticated queries to BioCyc
session = setup_biocyc_session(YOUR_MAIL, YOUR_PASSWORD)
# list of tuples (urls, gene, gene names, gpr, SGPR)
gprs = [getGPR(ec, session) for ec in ec_numbers]
```

### Session from environment variables

Email and password arguments can be specified outside of the script.

```bash
# this has to be the email that was used to create the account
export BIOCYC_EMAIL="YOUR EMAIL"
export BIOCYC_PASSWORD="YOUR STRONG PASSWORD"
```

And then, the script can be run as before:

```python
from gpr.auth_gpr import getGPR, setup_biocyc_session

ec_numbers = ["1.1.1.1", "1.1.1.204"]
# now without entering the email in password
#                             👇
session = setup_biocyc_session()
# list of tuples (urls, gene, gene names, gpr, SGPR)
gprs = [getGPR(ec, session) for ec in ec_numbers]
```
