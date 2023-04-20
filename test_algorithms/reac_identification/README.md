# Reac identification

This repository is a python module to evaulute the reaction identification algorithm (execute_jaccard)

```test_reac_identification.py``` is the main script for run the pytest and will use the following modules when run:
- Apply the reaction identification algorithm (functions_reac_identification.py.execute_jaccard)
- Other functions neccessary to run the test (functions_ast_gpr.py)

```reac_identification.py``` is the main script for evaluation and file result written.

## Usage

To run the pytest, use the following command:

```
pytest test_reac_identification.py
```

To write the results to file, use the following command:

```
python3 reac_identification.py
```
