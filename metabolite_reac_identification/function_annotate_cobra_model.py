from typing import Dict

import cobra
from cobra.io import write_sbml_model
import re

def annotate_cobra_model(model: cobra.Model, met_annotation: Dict, reac_annotation: Dict, out_file1: str =  "models/THG-beta1.1.1.xml", out_file2: str = "models/THG-beta1.1.xml"):
    """Update annotations of a model.
    Write two models to the file system.
    """
    # In[15]:
    with model:
        for metabolite in model.metabolites:
            if met_annotation.get(str(metabolite)[:-1]):
                a2 = met_annotation.get(str(metabolite)[:-1])
                for a in a2:
                    metabolite.annotation[a] = a2[a]
        for reaction, r5 in reac_annotation.items():
            model.reactions.get_by_id(reaction).annotation["kegg.reaction"] = r5[
                "kegg.reaction"
            ]
        write_sbml_model(model, out_file1)

    Output = open(out_file2, "w")
    tuple = (
        ('fbc:label="G_', 'fbc:label="'),
        ('name="G_', 'name="'),
        ('="meta_G_', '="'),
        ('="#?meta_[MGR]_', '="'),
        ('id="[MGR]_', 'id="'),
        ('species="M_', 'species="'),
        ('fbc:geneProduct="G_', 'fbc:geneProduct="'),
        ('idRef="R_', 'idRef="'),
        ('groups:name=".+?" ', ""),
        ('rdf:about="', 'rdf:about="#'),
    )

    with open(out_file1, "r") as model:
        model = model.read()
        for x, y in tuple:
            model = re.sub(x, y, model)
    Output.write(model)