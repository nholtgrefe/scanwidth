"""Configuration for edge-scanwidth reduction pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ReducerConfig:
    """Toggle reduction rules for :class:`scanwidth.edge_scanwidth.reduction.Reducer`.

    Parameters
    ----------
    use_sblocks : bool, optional
        If True, decompose into s-block subproblems.
    use_single_edge_rule : bool, optional
        If True, solve two-vertex single-edge blocks directly.
    use_single_root_cycle_rule : bool, optional
        If True, solve single-root cycle blocks directly.
    use_chain_suppression : bool, optional
        If True, apply chain-vertex suppression before delegating to solver.
    parallel_sblocks : bool, optional
        If True, solve independent s-block subproblems concurrently.
    sblock_max_workers : Optional[int], optional
        Maximum number of workers used for s-block parallelization.
        ``None`` delegates worker selection to the executor default.
    """

    use_sblocks: bool = True
    use_single_edge_rule: bool = True
    use_single_root_cycle_rule: bool = True
    use_chain_suppression: bool = True
    parallel_sblocks: bool = False
    sblock_max_workers: Optional[int] = None
