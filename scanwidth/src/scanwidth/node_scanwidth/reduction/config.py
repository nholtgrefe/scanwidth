"""Configuration for node-scanwidth reduction pipeline."""

from __future__ import annotations

from dataclasses import dataclass


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
    """

    use_sblocks: bool = True
    use_single_edge_shortcut: bool = True
    use_tree_shortcut: bool = True
    use_chain_suppression: bool = True
