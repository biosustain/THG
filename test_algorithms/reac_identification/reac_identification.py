import pandas as pd

reac = pd.read_csv('files/reac_results.csv')[['reac', 'expected_ec', 'message']]

x = reac[reac.message.isnull()].reset_index(drop=True)
x['ec'] = x.expected_ec
x['result'] = 'reaction correctly predicted'

y = reac[reac.message.str.contains("assert '", na=False)]
y = y.assign(ec=y.message.str.split('in (.+)',expand=True, regex=True)[1].str.replace("'",'',regex=False).str.replace("[",'', regex=False).str.replace("]",'', regex=False))
y['result'] = 'reaction incorrectly predicted'

z = reac[reac.message.str.contains('Jaccard result', na=False)]
z.insert(3, 'results', 'reaction not included in the test')

v = reac[reac.message.str.contains("not ec number in human model", na=False)]
v.insert(3, 'results', 'reaction not included in the test')

a_X = pd.concat([x, y, z, v])
a_X.rename(columns={'reac':'kegg_reac_id'}, inplace=True)
a_X = a_X.drop(columns=["message"])

a_X.to_excel('files/reac_identification.xlsx', engine='xlsxwriter', index=False)