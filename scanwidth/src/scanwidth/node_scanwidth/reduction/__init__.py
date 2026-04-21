"""Reduction helpers and orchestration for node scanwidth."""

from scanwidth.node_scanwidth.reduction.config import ReducerConfig
from scanwidth.node_scanwidth.reduction.reducer import Reducer

__all__ = [
    "Reducer",
    "ReducerConfig",
]
