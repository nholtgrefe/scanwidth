"""Public API entry point for node-scanwidth algorithms."""

from __future__ import annotations

from typing import Tuple

from scanwidth.dag import DAG
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.reduction.config import ReducerConfig
from scanwidth.node_scanwidth.reduction.reducer import Reducer
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth.node_scanwidth.solver.exact.exhaustive import ExhaustiveSolver
from scanwidth.node_scanwidth.solver.exact.ilp import ILPSolver
from scanwidth.node_scanwidth.solver.heuristic.greedy import GreedySolver
from scanwidth.node_scanwidth.solver.heuristic.random import RandomSolver
from scanwidth.node_scanwidth.solver.heuristic.simulated_annealing import (
    SimulatedAnnealingSolver,
)


def node_scanwidth(
    dag: DAG,
    algorithm: str = "exhaustive",
    reduce: bool = True,
    **kwargs: object,
) -> Tuple[int, Extension]:
    """Compute node scanwidth of a DAG using the selected algorithm."""
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
    """Instantiate a node-scanwidth solver from name and kwargs."""
    if algorithm == "exhaustive":
        return ExhaustiveSolver()
    if algorithm == "ilp":
        raw_time_limit = kwargs.pop("time_limit", None)
        return ILPSolver(
            time_limit=None if raw_time_limit is None else float(raw_time_limit),
            mip_rel_gap=float(kwargs.pop("mip_rel_gap", 0.0)),
            verbose=bool(kwargs.pop("verbose", False)),
        )
    if algorithm == "greedy":
        return GreedySolver()
    if algorithm == "random":
        return RandomSolver(seed=int(kwargs.pop("seed", 42)))
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
        "algorithm must be one of {'exhaustive', 'ilp', 'greedy', 'random', "
        "'simulated_annealing'}"
    )
