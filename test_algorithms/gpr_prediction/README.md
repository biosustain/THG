# GPR prediction

This repository is a python module to evalute the S-GPR building algorithm (getGPR)

```test_gpr_prediction.py``` is the main script for run the pytest and will use the following modules when run:
- Apply the S-GPR building algorithm (functions_auth_gpr.getGPR)
- Other functions neccessary to run the test (functions_ast_gpr.py)

```gpr_prediction.py``` is the main script for evaluation and file result written.

## Usage

To run the pytest, use the following command:

```
pytest test_gpr_prediction.py
```

To write the results to file, use the following command:

```
python3 gpr_prediction.py
```
