"""Metabolic task data model."""
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional

import cobra


class TaskResult(Enum):
    """Store whether a task was not run, passed, or failed."""

    NOT_RUN = auto()
    SUCCESS = auto()
    FAIL = auto()


@dataclass
class MetabolicTask:
    """Task extracted from TSV file.

    References
    ----------
    `Rasmus et al., 2014 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4017677/`_
    `RAVEN reference implementation https://github.com/SysBioChalmers/RAVEN/blob/800362e71b9d04ae2e2a16960b9d19aff346d37e/core/parseTaskList.m`_ .

    Parameters
    ----------
    id: str
    desc: str
    should_fail: bool
    uptaken: Optional[list[cobra.Reaction]]
        from IN field (since "in" is a keyword). Exchanges that must be uptaken.
        Note that this assumes canonical directionality of exchanges in the defaults.
    uptaken_lb: Optional[list[float]]
        Lower bound of uptakes. Default -1000.
    uptaken_ub: Optional[list[float]]
        Upper bound of uptakes. Default -1e4 (force directionality).
    out: Optional[list[cobra.Reaction]]
        from OUT field. Exchanges that must be producing.
    out_lb: Optional[list[float]]
        Default 1e-4 (force directionality).
    out_ub: Optional[list[float]]
        Default 1000.
    added_reaction: Optional[cobra.Reaction]
        from "EQU" field. If the reaction is present in the model, it is simply modified.
        The default of the bounds is unbounded reversibility (-1000, 1000).
    added_reaction_lb: float
        Default -1000.
    added_reaction_ub: float
        Default 1000.
    changed_reaction: Optional[str]
        Not implemented for now.
    changed_reaction_lb: float
    changed_reaction_ub: float
    comments: str
    references: Optional[str]
    explanations: str
    result: TaskResult
        Field to store the outcome of the task.
    allow_all_in: bool, default = False
        if ALLMETSIN[e] is passed as uptaken ("IN"), open all lower bound exchanges.
    allow_all_out: bool, default = False
        if ALLMETSIN[e] is passed as produced ("OUT"), open all upper bound exchanges.
    """

    id: str
    desc: str
    should_fail: bool
    uptaken: Optional[list[cobra.Reaction]]
    uptaken_lb: Optional[list[float]]
    uptaken_ub: Optional[list[float]]
    out: Optional[list[cobra.Reaction]]
    out_lb: Optional[list[float]]
    out_ub: Optional[list[float]]
    added_reaction: Optional[cobra.Reaction]
    added_reaction_lb: float
    added_reaction_ub: float
    changed_reaction: Optional[str]
    changed_reaction_lb: float
    changed_reaction_ub: float
    comments: str
    references: Optional[str]
    explanations: Optional[str]
    allow_all_in: bool = False
    allow_all_out: bool = False
    result: TaskResult = TaskResult.NOT_RUN
