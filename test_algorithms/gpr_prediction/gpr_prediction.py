import pandas as pd

x = pd.read_csv("files/ec-number-results.csv")[["ecnumber","expected_gprs","message"]]

y = x[x.message.isnull()].reset_index(drop=True)
y['predicted_gpr'] = y['expected_gprs']
y['result'] = 'expected gpr and predicted gpr were right'

z = x[x.message.str.contains("assert '", na=False)]
v = x[x.message.str.contains("None", na=False)]

w = pd.concat([z, v])
w['predicted_gpr'] = w.message.str.extract(" in '(.+)'\n")
w = w.drop(columns=["message"])

a_X =['2.4.1.17','1.11.1.7'] 
b_X =['2.4.1.224','2.7.4.6','6.2.1.5','6.4.1.4']
b_Y =['1.1.1.211','1.1.1.50','1.1.5.3','1.14.13.70','1.3.99.3','1.6.99.3','2.3.1.154','2.4.1.37','2.7.1.113', '2.7.1.2','2.7.1.48', '2.7.10.1','3.1.1.23','3.5.1.3','3.6.1.52', '3.6.3.14', '3.6.4.13', '4.2.1.3','4.2.3.2','6.2.1.17'] 
b_Z =['1.1.1.161','1.1.1.72','1.13.11.42','1.14.13.95','1.14.99.7','1.2.1.21','1.2.1.29','1.2.1.40','1.3.3.2','1.4.1.4','2.4.2.2','1.1.1.35','1.1.1.42','1.2.1.24','1.2.1.31','2.7.7.23','3.1.1.3','3.1.2.2','3.1.4.11','3.1.4.17','3.2.1.35','3.4.11.1','3.6.1.1','4.6.1.2']

a_X = z[z.ecnumber.isin(a_X)] # manually inspected to be right with the information available (incomplete Biocyc gene information)
b_X = z[z.ecnumber.isin(b_X)] # manually inspected to be wrong (the OR/AND operater)
b_Y = z[z.ecnumber.isin(b_Y)] # manually inspected to be wrong (the predicted gpr was incomplete)
b_Z = w[w.ecnumber.isin(b_Z)] # manually inspected to not to be None (predicted to be None)

ab_XXYZ = pd.concat([a_X,b_X,b_Y,b_Z])

a_Y = z[~z.ecnumber.isin(ab_XXYZ.ecnumber)] # manually inspected to be right 
a_Z = v[~v.ecnumber.isin(b_Z.ecnumber)] # # manually inspected to be right (None)

y_X = pd.concat([a_X.reset_index(drop=True), a_Y.reset_index(drop=True), a_Z.reset_index(drop=True)])
y_X ['result'] = 'predicted gpr was right'

y_Y = pd.concat([b_X.reset_index(drop=True), b_Y.reset_index(drop=True), b_Z.reset_index(drop=True)])
y_Y['result'] = 'expected and predicted gpr were wrong'

Y =pd.concat([y.reset_index(drop=True), y_X, y_Y])
Y = Y.drop(columns=["message"])
Y.rename(columns={'expected_gprs':'expected_gpr'}, inplace=True)

Y.to_excel('files/gpr_prediction.xlsx', engine='xlsxwriter', index=False)