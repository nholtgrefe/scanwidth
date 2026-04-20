"""Shared helper functions for node-scanwidth solvers."""

from __future__ import annotations

from typing import Set

import networkx as nx


def node_bag_size(graph: nx.DiGraph, connected_vertices: Set) -> int:
    """Return size of node-scanwidth bag for ``connected_vertices``."""
    return len(
        {
            u
            for (u, v) in graph.edges()
            if u not in connected_vertices and v in connected_vertices
        }
    )


def infinity_for(graph: nx.DiGraph) -> int:
    """Return a finite sentinel larger than any feasible scanwidth value."""
    return len(graph.nodes()) + len(graph.edges()) + 1
