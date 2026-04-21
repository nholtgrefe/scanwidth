"""Public API entry point for edge-scanwidth algorithms."""

from __future__ import annotations

from typing import Tuple

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.reduction.config import ReducerConfig
from scanwidth.edge_scanwidth.reduction.reducer import Reducer
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth.edge_scanwidth.solver.exact.exhaustive import ExhaustiveSolver
from scanwidth.edge_scanwidth.solver.exact.three_partition import (
    ThreePartitionSolver,
)
from scanwidth.edge_scanwidth.solver.exact.two_partition import TwoPartitionSolver
from scanwidth.edge_scanwidth.solver.exact.xp import XpSolver
from scanwidth.edge_scanwidth.solver.heuristic.cut_splitting import (
    CutSplittingSolver,
)
from scanwidth.edge_scanwidth.solver.heuristic.greedy import GreedySolver
from scanwidth.edge_scanwidth.solver.heuristic.random import RandomSolver
from scanwidth.edge_scanwidth.solver.heuristic.simulated_annealing import (
    SimulatedAnnealingSolver,
)
from scanwidth.extension import Extension


def edge_scanwidth(
    dag: DAG,
    algorithm: str = "xp",
    reduce: bool = True,
    **kwargs: object,
) -> Tuple[int, Extension]:
    """Compute edge scanwidth of a DAG using the selected algorithm.

    Parameters
    ----------
    dag : DAG
        Input directed acyclic graph.
    algorithm : str, optional
        Solver name. Supported values are:

        - ``"xp"``: exact XP algorithm with increasing ``k`` (default).
        - ``"exhaustive"``: exact exhaustive search over all extensions.
        - ``"two_partition"``: exact two-partition recursion.
        - ``"three_partition"``: exact 3-partition recursion.
        - ``"greedy"``: greedy heuristic.
        - ``"random"``: random extension.
        - ``"cut_splitting"``: recursive cut-splitting heuristic.
        - ``"simulated_annealing"``: simulated-annealing heuristic.

    reduce : bool, optional
        If True, apply s-block reduction before solving. Default is True.
    **kwargs : object
        Algorithm-specific keyword arguments forwarded to the matching
        solver class.

    Returns
    -------
    Tuple[int, Extension]
        Edge-scanwidth value and the corresponding extension.
    """
    solver = _build_solver(algorithm, kwargs)
    reducer_config = kwargs.pop("reducer_config", None)
    if kwargs:
        raise ValueError(
            f"Unexpected keyword arguments for algorithm '{algorithm}': "
            f"{sorted(kwargs)}"
        )
    if reducer_config is not None and not isinstance(reducer_config, ReducerConfig):
        raise TypeError("reducer_config must be a ReducerConfig instance.")

    if reduce:
        result = Reducer(
            config=ReducerConfig() if reducer_config is None else reducer_config,
        ).reduce_and_solve(dag=dag, solver=solver)
    else:
        result = solver.solve(dag)

    return result.value, result.extension


def _build_solver(algorithm: str, kwargs: dict) -> Solver:
    """Instantiate a solver from ``algorithm`` name and pop its kwargs."""
    if algorithm == "xp":
        if "k" in kwargs:
            raise ValueError(
                "Public edge_scanwidth(..., algorithm='xp') does not support "
                "fixed-k mode."
            )
        return XpSolver()
    if algorithm == "exhaustive":
        return ExhaustiveSolver()
    if algorithm == "two_partition":
        return TwoPartitionSolver()
    if algorithm == "three_partition":
        return ThreePartitionSolver()
    if algorithm == "greedy":
        return GreedySolver()
    if algorithm == "random":
        return RandomSolver(seed=int(kwargs.pop("seed", 42)))
    if algorithm == "cut_splitting":
        return CutSplittingSolver()
    if algorithm == "simulated_annealing":
        return SimulatedAnnealingSolver(
            max_iter=int(kwargs.pop("max_iter", 100)),
            p_in=float(kwargs.pop("p_in", 0.9)),
            p_stop=float(kwargs.pop("p_stop", 0.01)),
            init_ext=kwargs.pop("init_ext", "greedy"),  # type: ignore[arg-type]
            verbose=bool(kwargs.pop("verbose", True)),
            seed=int(kwargs.pop("seed", 42)),
        )
    raise ValueError(
        "algorithm must be one of {'xp', 'exhaustive', 'two_partition', "
        "'greedy', 'random', 'cut_splitting', 'simulated_annealing'}"
    )
