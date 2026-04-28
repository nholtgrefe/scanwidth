"""Exact XP solver (increasing-k) for node scanwidth."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Union

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth._utils import infinity_for


class RpswTable:
    """Memoization table for restricted partial scanwidth (rpsw)."""

    def __init__(self) -> None:
        self._entries: Dict[str, Tuple[int, List]] = {}

    def __contains__(self, key: str) -> bool:
        return key in self._entries

    def __getitem__(self, key: str) -> Tuple[int, List]:
        return self._entries[key]

    def __setitem__(self, key: str, value: Tuple[int, List]) -> None:
        self._entries[key] = value

    def get(self, key: str) -> Optional[Tuple[int, List]]:
        """Return the cached value for ``key`` or ``None``."""
        return self._entries.get(key)

    def drop_infeasible(self, infinity: int) -> None:
        """Remove entries whose cached rpsw equals ``infinity``."""
        for key in list(self._entries):
            if self._entries[key][0] == infinity:
                del self._entries[key]


def _delta_in_distinct_parents(
    graph: nx.DiGraph,
    vertex_set: Union[Set, List],
) -> int:
    """Return number of distinct parents outside ``vertex_set`` with an edge into it."""
    vertices = set(vertex_set)
    boundary_parents = {
        u for (u, v) in graph.edges() if v in vertices and u not in vertices
    }
    return len(boundary_parents)


@dataclass(frozen=True)
class XpSolver(Solver):
    """Exact XP node-scanwidth solver using increasing ``k``."""

    k: Optional[int] = None

    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth exactly using the XP algorithm."""
        graph = dag.graph
        infinity = infinity_for(graph)
        table = RpswTable()

        if self.k is not None:
            return self._solve_fixed_k(graph, self.k, table, infinity)

        for k in range(1, infinity + 1):
            sw, sigma = self._restricted_partial_scanwidth(
                graph, graph.nodes(), table, infinity, k,
            )
            if sw == infinity:
                table.drop_infeasible(infinity)
                continue
            return SolverResult(value=sw, extension=Extension(dag, sigma))

        raise RuntimeError("XP solver failed to converge.")

    def _solve_fixed_k(
        self,
        graph: nx.DiGraph,
        k: int,
        table: RpswTable,
        infinity: int,
    ) -> SolverResult:
        """Solve the fixed-``k`` variant, raising if no solution exists."""
        sw, sigma = self._restricted_partial_scanwidth(
            graph, graph.nodes(), table, infinity, k,
        )
        if sw == infinity:
            raise ValueError(
                f"No extension exists with node scanwidth <= {k}.",
            )
        return SolverResult(value=sw, extension=Extension(DAG(graph), sigma))

    def _restricted_partial_scanwidth(
        self,
        graph: nx.DiGraph,
        vertices: Union[Set, List],
        table: RpswTable,
        infinity: int,
        k: int,
    ) -> Tuple[int, List]:
        """Compute restricted partial scanwidth with component splitting."""
        vertex_list = list(vertices)
        subgraph = graph.subgraph(vertex_list)
        roots = [v for v in subgraph.nodes() if subgraph.in_degree(v) == 0]

        key = repr(sorted(roots))
        cached = table.get(key)
        if cached is not None:
            return cached

        delta_in_W = _delta_in_distinct_parents(graph, vertices)
        rpsw = infinity
        sigma: List = []

        if len(vertex_list) == 1 and delta_in_W <= k:
            rpsw = delta_in_W
            sigma = [vertex_list[0]]
        else:
            components = list(nx.weakly_connected_components(subgraph))

            if len(components) > 1:
                rpsw_list = []
                for U_i in components:
                    rpsw_i, sigma_i = self._restricted_partial_scanwidth(
                        graph, U_i, table, infinity, k,
                    )
                    rpsw_list.append(rpsw_i)
                    sigma = sigma + sigma_i
                rpsw = max(rpsw_list)
            elif len(vertex_list) > 1 and delta_in_W <= k:
                for rho in roots:
                    new_vertices = vertex_list.copy()
                    new_vertices.remove(rho)
                    rpsw1, sigma1 = self._restricted_partial_scanwidth(
                        graph, new_vertices, table, infinity, k,
                    )
                    rpsw_prime = max(rpsw1, delta_in_W)
                    if rpsw_prime < rpsw:
                        rpsw = rpsw_prime
                        sigma = sigma1 + [rho]

        table[key] = (rpsw, sigma)
        return rpsw, sigma
