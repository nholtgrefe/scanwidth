"""Solver classes for edge scanwidth."""

from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth.edge_scanwidth.solver.exact.exhaustive import BruteForceSolver
from scanwidth.edge_scanwidth.solver.exact.three_partition import (
    ThreePartitionSolver,
)
from scanwidth.edge_scanwidth.solver.exact.two_partition import TwoPartitionSolver
from scanwidth.edge_scanwidth.solver.exact.xp import RpswTable, XpSolver
from scanwidth.edge_scanwidth.solver.heuristic.cut_splitting import (
    CutSplittingSolver,
)
from scanwidth.edge_scanwidth.solver.heuristic.greedy import GreedySolver
from scanwidth.edge_scanwidth.solver.heuristic.random import RandomSolver
from scanwidth.edge_scanwidth.solver.heuristic.simulated_annealing import (
    SimulatedAnnealingSolver,
)

__all__ = [
    "Solver",
    "XpSolver",
    "RpswTable",
    "BruteForceSolver",
    "TwoPartitionSolver",
    "ThreePartitionSolver",
    "GreedySolver",
    "RandomSolver",
    "CutSplittingSolver",
    "SimulatedAnnealingSolver",
]
