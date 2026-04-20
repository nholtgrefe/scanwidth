"""Base solver interfaces for edge-scanwidth algorithms."""

from __future__ import annotations

from abc import ABC, abstractmethod

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult


class Solver(ABC):
    """Abstract base class for edge-scanwidth solvers."""

    @abstractmethod
    def solve(self, dag: DAG) -> SolverResult:
        """Solve edge scanwidth on the provided DAG.

        Parameters
        ----------
        dag : DAG
            Input graph instance.

        Returns
        -------
        SolverResult
            Standardized solver result.
        """
        raise NotImplementedError
