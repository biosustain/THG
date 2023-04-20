# GPR prediction

This repository is a python module to evalute the S-GPR building algorithm (functions_auth_gpr.getGPR)

"test_gpr_prediction.py" is the script to execute the getGPR function and will use the following modules when run:
- functions_auth_gpr
- functions_ast_gpr

The result file can be generated from "gpr_prediction.py"

## Usage

The test is run using pytest:

```
pytest test_gpr_prediction.py
```

To generate the result file, use the following command:

```
python3 gpr_prediction.py
```
