"""Shared types for edge-scanwidth solvers and APIs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from scanwidth.extension import Extension


@dataclass(frozen=True)
class SolverResult:
    """Container with standardized solver output.

    Parameters
    ----------
    value : int
        Computed edge-scanwidth value.
    extension : Extension
        Extension corresponding to the returned scanwidth.
    history : Optional[List[int]], optional
        Optional optimization trajectory (for example for annealing).
    diagnostics : Dict[str, Any], optional
        Extra metadata produced by the solver.
    """

    value: int
    extension: Extension
    history: Optional[List[int]] = None
    diagnostics: Dict[str, Any] = field(default_factory=dict)
