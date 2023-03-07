"""Test metabolic tasks.

Memote custom tests are defined at the beginning,
followed by the logic of the tests.
"""

from os.path import dirname, join

import cobra
import pandas as pd
import pytest
from memote.utils import annotate, truncate, wrapper
from metabolic_tasks import iterate_over_tasks, TaskResult


@pytest.mark.parametrize(
    "task_group", ["Essential", "CellfieConsensus", "Full", "VerifyModel"]
)
@annotate(
    title="Metabolic tasks",
    format_type="percent",
    message=dict(),
    data=dict(),
    metric=dict(),
)
def test_essential_metabolic_tasks(model: cobra.Model, task_group: str):
    """Count number of passing metabolic tasks."""
    task_file = f"metabolicTasks_{task_group}.txt"
    task_input = pd.read_csv(
        join(dirname(__file__), "data", task_file),
        sep="\t",
        comment="#",
    )
    tasks = iterate_over_tasks(model, task_input)
    ann = test_essential_metabolic_tasks.annotation
    ann["data"][task_group] = [
        f"{task.desc} ({task.id})" for task in tasks if task.result == TaskResult.FAIL
    ]
    ann["metric"][task_group] = len(
        [task for task in tasks if task.result == TaskResult.SUCCESS]
    ) / len(tasks)
    ann["message"][task_group] = wrapper.fill(
        f"""Passed {task_group} metabolic tasks.

        {ann["metric"][task_group]:.2%} of the {len(tasks)} tasks where succesful. Failed tasks: {truncate(ann["data"][task_group])}"""
    )
    assert ann["metric"][task_group] == 1.0, ann["message"][task_group]


if __name__ == "__main__":
    df = pd.read_csv(
        join(dirname(__file__), "data", "metabolicTasks_VerifyModel.txt"),
        sep="\t",
        comment="#",
    )
    model = cobra.io.read_sbml_model(join(dirname(__file__), "..", "models/THG-2023-02-25.xml")) # other models in models folder:Human-GEM_2022-06-21.xml, THG-beta1.xml, THG-beta2.xml
    tasks = iterate_over_tasks(model, df)
    print(f"Failed -> {[task.id for task in tasks if task.result == TaskResult.FAIL]}")
    print(
        f"Success -> {[task.id for task in tasks if task.result == TaskResult.SUCCESS]}"
    )
