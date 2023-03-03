from functions_compare_models import *


model = read_sbml_model('models/THG-2023-02-25.xml') # new model
model2 = read_sbml_model('models/Human-GEM_2022-06-21.xml') # org. model


Output = xlsxwriter.Workbook("reports/THG_vs_Human1.xlsx")


MODEL = Output.add_worksheet("MODEL")
GROUPS = Output.add_worksheet("GROUPS")
COMPS = Output.add_worksheet("COMPS")
METS = Output.add_worksheet("METS")
METS2 = Output.add_worksheet("METS2")
RXNS = Output.add_worksheet("RXNS")
GENES = Output.add_worksheet("GENES")


defmodel(model, MODEL, Output, model2)
defgroup(GROUPS,Output)
compartment(model, COMPS, Output, model2)
#metabolite(model, METS, Output, model2)
metabolite_2(model, METS, Output, model2) # Only for merged models because it accounts for lipidbank annotation
#metabolite2(model, model2, METS2,Output)
metabolite2_2(model, model2, METS2,Output) # Only for merged models because it accounts for lipidbank annotation 
reaction(model, RXNS, Output, model2)
gene(model2, GENES, Output, model)


Output.close()
 
