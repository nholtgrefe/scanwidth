"""Configuration for node-scanwidth reduction pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ReducerConfig:
    """Toggle reduction rules for :class:`scanwidth.node_scanwidth.reduction.Reducer`.

    Parameters
    ----------
    use_sblocks : bool, optional
        If True, decompose into s-block subproblems.
    use_single_edge_shortcut : bool, optional
        If True, solve two-vertex single-edge blocks directly.
    use_tree_shortcut : bool, optional
        If True, solve directed rooted trees directly.
    use_chain_suppression : bool, optional
        If True, apply chain-parent suppression on blocks before solving.
    use_reticulation_path_shortcut : bool, optional
        If True, apply the reticulation-path closed-form shortcut on
        eligible non-root s-blocks after chain suppression.
    parallel_sblocks : bool, optional
        If True, solve independent s-block subproblems concurrently.
    sblock_max_workers : Optional[int], optional
        Maximum number of workers used for s-block parallelization.
        ``None`` delegates worker selection to the executor default.
    """

    use_sblocks: bool = True
    use_single_edge_shortcut: bool = True
    use_tree_shortcut: bool = True
    use_chain_suppression: bool = True
    use_reticulation_path_shortcut: bool = True
    parallel_sblocks: bool = False
    sblock_max_workers: Optional[int] = None
