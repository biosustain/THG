import pandas as pd


def process_met_file(met_file):
    
    m = pd.read_csv(met_file)[["identifier", "name", "formula", "expected", "message"]]
    m["identifier"] = m["identifier"].replace("[a-z]+","", regex=True)

    x = m[m["message"].str.contains("assert '", na=False)]
    x.insert(5, 'results', 'metabolite incorrectly annotated')

    y = m[m['message'].isnull()].reset_index(drop=True)
    y['result'] = 'metabolite correctly annotated'

    z = m[m["message"].str.contains("None", na=False)]
    z.insert(5, 'results', 'metabolite not identified')

    x = pd.concat([x, y])
    x['id'] = x['message'].str.split('\\\\t',expand=True, regex=True)[6]
    x['id'] = x['id'].fillna(x['expected'])
    x = x.drop(columns=["message"])

    a_x = pd.concat([x, z])
    a_x = a_x.drop(columns=["message"]).reset_index(drop=True)
    a_x = a_x[['identifier','name','formula','expected','id','result']]
    a_x.rename(columns={'expected':'expected_id'}, inplace=True)
    
    return a_x
    
process_met_file('files/met_results_010.csv')

with pd.ExcelWriter('files/metabolite_identification.xlsx') as writer:
    process_met_file('files/met_results_010.csv').to_excel(writer, sheet_name='met_threshold_010', index=False)
    process_met_file('files/met_results_060.csv').to_excel(writer, sheet_name='met_threshold_060', index=False)
    process_met_file('files/met_results_070.csv').to_excel(writer, sheet_name='met_threshold_070', index=False)
    process_met_file('files/met_results_080.csv').to_excel(writer, sheet_name='met_threshold_080', index=False)
    process_met_file('files/met_results_082.csv').to_excel(writer, sheet_name='met_threshold_082', index=False)