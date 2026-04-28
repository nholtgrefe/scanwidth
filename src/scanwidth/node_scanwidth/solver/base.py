"""Base class for node-scanwidth solvers."""

from __future__ import annotations

from abc import ABC, abstractmethod

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult


class Solver(ABC):
    """Abstract interface for node-scanwidth solver implementations."""

    @abstractmethod
    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth on ``dag`` and return a standardized result."""
