"""Exact XP solver (increasing-k) for node scanwidth."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Dict, Hashable, List, Optional, Set, Tuple, Union

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth._utils import (
    delta_in_parents,
    induced_subgraph_roots,
    induced_weakly_connected_components,
    infinity_for,
)


class RpswTable:
    """Memoization table for restricted partial scanwidth (rpsw)."""

    def __init__(self) -> None:
        self._entries: Dict[Hashable, Tuple[int, List]] = {}

    def __contains__(self, key: Hashable) -> bool:
        return key in self._entries

    def __getitem__(self, key: Hashable) -> Tuple[int, List]:
        return self._entries[key]

    def __setitem__(self, key: Hashable, value: Tuple[int, List]) -> None:
        self._entries[key] = value

    def get(self, key: Hashable) -> Optional[Tuple[int, List]]:
        """Return the cached value for ``key`` or ``None``."""
        return self._entries.get(key)

    def drop_infeasible(self, infinity: int) -> None:
        """Remove entries whose cached rpsw equals ``infinity``."""
        for key in list(self._entries):
            if self._entries[key][0] == infinity:
                del self._entries[key]


@dataclass(frozen=True)
class XpSolver(Solver):
    """Exact XP node-scanwidth solver using increasing ``k``."""

    k: Optional[int] = None

    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth exactly using the XP algorithm."""
        graph = dag.graph
        aux_graph, represented_vertices, expansion_history = (
            self._suppress_auxiliary_chain_vertices(graph)
        )
        effective_repr: Optional[Dict[object, Set]]
        if expansion_history:
            effective_repr = represented_vertices
        else:
            # No suppression happened: avoid set-expansion overhead in recursion.
            effective_repr = None
        infinity = infinity_for(aux_graph)
        table = RpswTable()

        if self.k is not None:
            return self._solve_fixed_k(
                original_graph=graph,
                aux_graph=aux_graph,
                represented_vertices=effective_repr,
                expansion_history=expansion_history,
                k=self.k,
                table=table,
                infinity=infinity,
            )

        for k in range(1, infinity + 1):
            sw, sigma = self._restricted_partial_scanwidth(
                original_graph=graph,
                aux_graph=aux_graph,
                vertices=aux_graph.nodes(),
                represented_vertices=effective_repr,
                table=table,
                infinity=infinity,
                k=k,
            )
            if sw == infinity:
                table.drop_infeasible(infinity)
                continue
            sigma = self._expand_auxiliary_ordering(sigma, expansion_history)
            return SolverResult(value=sw, extension=Extension(dag, sigma))

        raise RuntimeError("XP solver failed to converge.")

    def _solve_fixed_k(
        self,
        original_graph: nx.DiGraph,
        aux_graph: nx.DiGraph,
        represented_vertices: Optional[Dict[object, Set]],
        expansion_history: List[Tuple[object, object]],
        k: int,
        table: RpswTable,
        infinity: int,
    ) -> SolverResult:
        """Solve the fixed-``k`` variant, raising if no solution exists."""
        sw, sigma = self._restricted_partial_scanwidth(
            original_graph=original_graph,
            aux_graph=aux_graph,
            vertices=aux_graph.nodes(),
            represented_vertices=represented_vertices,
            table=table,
            infinity=infinity,
            k=k,
        )
        if sw == infinity:
            raise ValueError(
                f"No extension exists with node scanwidth <= {k}.",
            )
        sigma = self._expand_auxiliary_ordering(sigma, expansion_history)
        return SolverResult(
            value=sw,
            extension=Extension(DAG(original_graph), sigma),
        )

    def _restricted_partial_scanwidth(
        self,
        original_graph: nx.DiGraph,
        aux_graph: nx.DiGraph,
        vertices: Union[Set, List],
        represented_vertices: Optional[Dict[object, Set]],
        table: RpswTable,
        infinity: int,
        k: int,
    ) -> Tuple[int, List]:
        """Compute restricted partial scanwidth with component splitting."""
        vertex_list = list(vertices)
        roots = induced_subgraph_roots(aux_graph, vertex_list)

        key = tuple(sorted(roots, key=str))
        cached = table.get(key)
        if cached is not None:
            return cached

        if represented_vertices is None:
            delta_in_W = delta_in_parents(aux_graph, vertex_list, sink=True)
        else:
            expanded_vertex_set = self._expanded_vertex_set(
                vertex_list,
                represented_vertices,
            )
            delta_in_W = delta_in_parents(
                original_graph,
                expanded_vertex_set,
                sink=True,
            )
        rpsw = infinity
        sigma: List = []

        if len(vertex_list) == 1 and delta_in_W <= k:
            rpsw = delta_in_W
            sigma = [vertex_list[0]]
        else:
            components = induced_weakly_connected_components(aux_graph, vertex_list)

            if len(components) > 1:
                rpsw_list = []
                for U_i in components:
                    rpsw_i, sigma_i = self._restricted_partial_scanwidth(
                        original_graph=original_graph,
                        aux_graph=aux_graph,
                        vertices=U_i,
                        represented_vertices=represented_vertices,
                        table=table,
                        infinity=infinity,
                        k=k,
                    )
                    rpsw_list.append(rpsw_i)
                    sigma.extend(sigma_i)
                rpsw = max(rpsw_list)
            elif len(vertex_list) > 1 and delta_in_W <= k:
                for rho in roots:
                    new_vertices = vertex_list.copy()
                    new_vertices.remove(rho)
                    rpsw1, sigma1 = self._restricted_partial_scanwidth(
                        original_graph=original_graph,
                        aux_graph=aux_graph,
                        vertices=new_vertices,
                        represented_vertices=represented_vertices,
                        table=table,
                        infinity=infinity,
                        k=k,
                    )
                    rpsw_prime = max(rpsw1, delta_in_W)
                    if rpsw_prime < rpsw:
                        rpsw = rpsw_prime
                        sigma = sigma1 + [rho]

        table[key] = (rpsw, sigma)
        return rpsw, sigma

    @staticmethod
    def _expanded_vertex_set(
        aux_vertices: List,
        represented_vertices: Dict[object, Set],
    ) -> Set:
        """Expand an auxiliary sink-set to original graph vertices.

        Parameters
        ----------
        aux_vertices : List
            Vertices in the auxiliary graph.
        represented_vertices : Dict[object, Set]
            Map from an auxiliary vertex to original represented vertices.

        Returns
        -------
        Set
            Expanded sink-set in original-graph coordinates.
        """
        expanded: Set = set()
        for vertex in aux_vertices:
            expanded.update(represented_vertices[vertex])
        return expanded

    @staticmethod
    def _suppress_auxiliary_chain_vertices(
        graph: nx.DiGraph,
    ) -> Tuple[nx.DiGraph, Dict[object, Set], List[Tuple[object, object]]]:
        """Suppress all degree-(1,1) vertices in an auxiliary graph.

        Parameters
        ----------
        graph : nx.DiGraph
            Original graph.

        Returns
        -------
        Tuple[nx.DiGraph, Dict[object, Set], List[Tuple[object, object]]]
            Auxiliary graph, represented-vertex map, and expansion history.
        """
        aux_graph = graph.copy()
        represented_vertices: Dict[object, Set] = {v: {v} for v in aux_graph.nodes()}
        expansion_history: List[Tuple[object, object]] = []
        candidates = deque(
            vertex
            for vertex in aux_graph.nodes()
            if aux_graph.in_degree(vertex) == 1 and aux_graph.out_degree(vertex) == 1
        )
        while candidates:
            vertex = candidates.popleft()
            if vertex not in aux_graph:
                continue
            if aux_graph.in_degree(vertex) != 1 or aux_graph.out_degree(vertex) != 1:
                continue
            predecessor = next(iter(aux_graph.predecessors(vertex)))
            successor = next(iter(aux_graph.successors(vertex)))
            aux_graph.remove_node(vertex)
            if predecessor != successor:
                aux_graph.add_edge(predecessor, successor)
            represented_vertices[successor].update(represented_vertices[vertex])
            del represented_vertices[vertex]
            expansion_history.append((successor, vertex))

            # Local degree changes can create new suppressible vertices.
            candidates.append(predecessor)
            candidates.append(successor)

        return aux_graph, represented_vertices, expansion_history

    @staticmethod
    def _expand_auxiliary_ordering(
        ordering: List,
        expansion_history: List[Tuple[object, object]],
    ) -> List:
        """Expand suppressed vertices back into an extension ordering.

        Parameters
        ----------
        ordering : List
            Ordering over auxiliary graph vertices.
        expansion_history : List[Tuple[object, object]]
            Pairs ``(kept_vertex, suppressed_vertex)`` in suppression order.

        Returns
        -------
        List
            Ordering over original graph vertices.
        """
        if not expansion_history:
            return list(ordering)

        next_vertex: Dict[object, Optional[object]] = {}
        for idx, vertex in enumerate(ordering):
            next_vertex[vertex] = ordering[idx + 1] if idx + 1 < len(ordering) else None

        for kept_vertex, suppressed_vertex in reversed(expansion_history):
            next_vertex[suppressed_vertex] = next_vertex[kept_vertex]
            next_vertex[kept_vertex] = suppressed_vertex

        if not ordering:
            return []
        restored: List = []
        current: Optional[object] = ordering[0]
        while current is not None:
            restored.append(current)
            current = next_vertex[current]
        return restored
