# Algorithm tests

This repository is a python module for the algorithm validation.

The tests and the results of the algorithm validation can be found in the different folders.

1. gpr_prediction: Assessing the S-GPR building algorithm (getGPR).

2. mass_balance: Assessing the algorithm to mass balance metabolic reactions (mass_balance).

3. metabolite_identification: Assessing the text similarity algorithm to identify metabolites (identify_metabolite).

4. reac_identification: Assessing the JI-based algorithm to identify metabolic reactions (execute_jaccard).

Each folder contains
- the functions (or equations) neccessary to perform the tests.
- the tests.
- the code to assess the algorithm test's results (folder_name.py).
- files (including result file (folder_name.xlsx)) and models in the respective subfolders.

To run each test:

```bash
pytest test_folder_name.py
```
