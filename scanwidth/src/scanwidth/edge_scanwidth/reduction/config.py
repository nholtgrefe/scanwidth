"""Configuration for edge-scanwidth reduction pipeline."""

from __future__ import annotations

from dataclasses import dataclass


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
    use_flow_suppression : bool, optional
        If True, apply flow-node suppression before delegating to solver.
    """

    use_sblocks: bool = True
    use_single_edge_rule: bool = True
    use_single_root_cycle_rule: bool = True
    use_flow_suppression: bool = True
