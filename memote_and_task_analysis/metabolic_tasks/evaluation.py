"""Evaluation of Metabolic Tasks."""
import re

import cobra
import pandas as pd

from task import MetabolicTask, TaskResult

COMP_REGEX = re.compile(r"\[([a-z]+)\]")
ARROW_REGEX = re.compile(r"<=>|=>|<=|<->|<-|->")


def gather_mets_from_reaction_str(
    model: cobra.Model, reaction_str: str
) -> dict[cobra.Metabolite, int]:
    reactants_str, products_str = ARROW_REGEX.split(reaction_str)

    def gather_summation(half_reaction) -> dict[cobra.Metabolite, int]:
        # gather comparments
        mets_comp = [
            (COMP_REGEX.sub("", parsed_met).strip(), COMP_REGEX.search(parsed_met)[1])
            for parsed_met in half_reaction.split(" + ")
        ]
        # gather possible stoichimetry
        stoich = [int(met.split(" ")[0]) if " " in met else 1 for met, _ in mets_comp]
        mets_comp = [
            (met.split(" ")[1] if " " in met else met, comp) for met, comp in mets_comp
        ]

        return {
            next(
                met
                for met in model.metabolites
                if met.name == parsed_met
                # TODO: this is fallible
                and met.compartment == comp
            ): st
            for st, (parsed_met, comp) in zip(stoich, mets_comp)
        }

    try:
        reactants = gather_summation(reactants_str)
        products = gather_summation(products_str)
    except StopIteration:
        raise StopIteration(
            f"""No metabolite with names
             '{[COMP_REGEX.sub('', reactant).strip() for reactant in reactants_str.split(' + ')]}' -> '{[COMP_REGEX.sub('', product).strip() for product in products_str.split(' + ')]}'
             and comp '{[COMP_REGEX.search(reactant)[1] for reactant in reactants_str.split(' + ')]}' -> '{[COMP_REGEX.search(product)[1] for product in products_str.split(' + ')]}' """
        )
    metabolites = {met: -1 * st for met, st in reactants.items()}
    metabolites.update({met: st for met, st in products.items()})
    return metabolites


def read_metabolic_task(
    model: cobra.Model, df: pd.DataFrame
) -> tuple[int, MetabolicTask]:
    """Read a metabolic task, iterating over a :code:`pd.DataFrame`."""
    df_iterator = df.itertuples()
    row = next(df_iterator)
    # true for first iteration; false if any other iteration contains an ID
    iterating = True
    task = None
    i = 0
    uptaken, uptaken_lb, uptaken_ub = [], [], []
    produced, produced_lb, produced_ub = [], [], []
    all_mets_in = False
    all_mets_out = False
    while iterating:
        i += 1
        # IN metabolites, uptaken
        if pd.isna(row.IN):
            pass
        elif "ALLMETSIN[e]" in row.IN:
            all_mets_in = True
        else:
            uptaken_mets = row.IN.split(";")
            uptaken = [
                ex
                for uptaken_met in uptaken_mets
                for ex in model.exchanges
                if any(met.name == uptaken_met[:-3] for met in ex.metabolites)
            ]
            # switch since they assume opposite directionality
            lb = -1000 if pd.isna(row.IN_UB) else -row.IN_UB
            ub = -1e-04 if pd.isna(row.IN_LB) else -row.IN_LB
            uptaken_lb = [lb for _ in range(len(uptaken))]
            uptaken_ub = [ub for _ in range(len(uptaken))]
        # OUT metabolites, secreted
        if pd.isna(row.OUT):
            pass
        elif "ALLMETSIN[e]" in row.OUT:
            all_mets_out = True
        else:
            produced_mets = row.OUT.split(";")
            produced = [
                ex
                for produced_met in produced_mets
                for ex in model.exchanges
                if any(met.name == produced_met[:-3] for met in ex.metabolites)
            ]
            lb = 1e-04 if pd.isna(row.OUT_LB) else row.OUT_LB
            ub = 1000 if pd.isna(row.OUT_UB) else row.OUT_UB
            produced_lb = [lb for _ in range(len(produced))]
            produced_ub = [ub for _ in range(len(produced))]
        # Added reaction
        if pd.isna(row.EQU):
            added_reaction = None
        else:
            if row.EQU in model.reactions:
                added_reaction = model.reactions.get_by_id(row.EQU)
            else:
                added_reaction = cobra.Reaction(f"MT_ADDED_{i}")
                model.add_reactions([added_reaction])
                metabolites = gather_mets_from_reaction_str(model, row.EQU)
                added_reaction.add_metabolites(metabolites)
        added_reaction_lb, added_reaction_ub = (
            -1000 if pd.isna(row.EQU_LB) else row.EQU_LB,
            1000 if pd.isna(row.EQU_UB) else row.EQU_UB,
        )
        if not pd.isna(row.CHANGED_RXN):
            raise NotImplementedError("Changed reactions are not implemented.")
        if task is None:
            task = MetabolicTask(
                row.ID,
                row.DESCRIPTION,
                False if pd.isna(row.SHOULD_FAIL) else row.SHOULD_FAIL,
                uptaken,
                uptaken_lb,
                uptaken_ub,
                produced,
                produced_lb,
                produced_ub,
                added_reaction,
                added_reaction_lb,
                added_reaction_ub,
                None,
                0,
                0,
                row.COMMENTS,
                row.REFERENCES if "REFERENCES" in df.columns else None,
                row.EXPLANATIONS
                if "EXPLANATIONS" in df.columns
                else row.COMMENTS
                if "COMMENTS" in df.columns
                else None,
                all_mets_in,
                all_mets_out,
            )
        else:
            task.uptaken += uptaken
            task.uptaken_lb += uptaken_lb
            task.uptaken_ub += uptaken_ub
            task.out += produced
            task.out_lb += produced_lb
            task.out_ub += produced_ub
        try:
            row = next(df_iterator)
        except StopIteration:
            break
        iterating = pd.isna(row.ID)
    return i, task


def apply_metabolic_task(model: cobra.Model, task: MetabolicTask):
    """Apply a metabolic task to a model **in-place**."""
    if task.uptaken is not None:
        for exchange, lb, ub in zip(task.uptaken, task.uptaken_lb, task.uptaken_ub):
            model.reactions.get_by_id(exchange.id).bounds = lb, ub
    if task.out is not None:
        for exchange, lb, ub in zip(task.out, task.out_lb, task.out_ub):
            model.reactions.get_by_id(exchange.id).bounds = lb, ub
    if task.allow_all_in:
        for exchange in model.exchanges:
            exchange.lower_bound = -1000
    if task.allow_all_out:
        for exchange in model.exchanges:
            exchange.upper_bound = 1000

    if task.added_reaction is not None:
        reac_id = task.added_reaction.id
        lb, ub = task.added_reaction_lb, task.added_reaction_ub
        if reac_id in model.reactions:
            model.reactions.get_by_id(reac_id).bounds = lb, ub
        else:
            raise NotImplementedError(
                "Not handling reactions not in the model for now."
            )
    if task.changed_reaction is not None:
        model.reactions.get_by_id(task.changed_reaction).bounds = (
            task.changed_reaction_lb,
            task.changed_reaction_ub,
        )


def iterate_over_tasks(model: cobra.Model, df: pd.DataFrame) -> list[MetabolicTask]:
    """Extract and apply each task separately, gathering their results."""
    nrows = df.shape[0]
    tasks: list[MetabolicTask] = []
    df = df.copy()
    # change spaces so that they are valid identifiers and can be accessed from
    # the NamedTuple returned by df.itertuples
    df.columns = df.columns.str.replace(" ", "_")

    i = 0
    while i < nrows:
        with model as m:
            advanced, task = read_metabolic_task(m, df.iloc[i:, :])
            apply_metabolic_task(m, task)
            obj_value = m.slim_optimize()
            task.result = (
                TaskResult.FAIL
                if (obj_value is None or obj_value < 1e-4) and not task.should_fail
                else TaskResult.SUCCESS
            )
            tasks.append(task)
        i += advanced
    return tasks
