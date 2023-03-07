from function_metabolite_identification import *
from function_reac_identification import *
from function_annotate_cobra_model import *
import cobra.io
import numpy as np
import re


model = 'models/Human-GEM_2022-06-21.xml'
database = 'models/Human Database.xml'

cobra_model = cobra.io.read_sbml_model(model)


# Metabolites

metabolites = gather_metabolites(cobra_model)
generate_met_annotation(metabolites)
met_annotation = process_annotation()


# Reactions

reac = process_reac(model, '[A-Z]+[0-9]+[a-z]+[0-9]*', 'MAM02040', 'MAM02039')
reac = replace_met_id_by_met_kegg(reac, gather_kegg_metabolites(model))
reac_y = process_reac(database, '([A-Z][0-9]+_?[a-z]+[0-9]*)', 'C00080', 'C00001')
jaccard = execute_jaccard(reac_y, reac)
reac_annotation = process_jaccard(reac_y, reac, jaccard)


# SBML

annotate_cobra_model(cobra_model, met_annotation, reac_annotation)
