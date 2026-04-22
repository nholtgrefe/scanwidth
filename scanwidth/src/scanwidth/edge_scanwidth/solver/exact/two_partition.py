"""Two-partition exact solver for edge scanwidth."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set, Tuple, Union

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth._utils import delta_in, infinity_for
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension


@dataclass(frozen=True)
class TwoPartitionSolver(Solver):
    """Exact solver using recursive two-partition decomposition.

    This corresponds to the component-splitting restricted-partial-scanwidth
    recursion without memoization.
    """

    def solve(self, dag: DAG) -> SolverResult:
        """Solve edge scanwidth exactly using two-partition recursion.

        Parameters
        ----------
        dag : DAG
            Input graph instance.

        Returns
        -------
        SolverResult
            Standardized solver result.
        """
        graph = dag.graph
        infinity = infinity_for(graph)
        sw, sigma = self._restricted_partial_scanwidth(
            graph=graph,
            vertices=set(graph.nodes()),
            k=infinity,
            infinity=infinity,
        )
        return SolverResult(value=sw, extension=Extension(dag, sigma))

    def _restricted_partial_scanwidth(
        self,
        graph: nx.DiGraph,
        vertices: Union[Set, List],
        k: int,
        infinity: int,
    ) -> Tuple[int, List]:
        """Compute restricted partial scanwidth with component splitting."""
        vertex_list = list(vertices)
        subgraph = graph.subgraph(vertex_list)
        roots = [v for v in subgraph.nodes() if subgraph.in_degree(v) == 0]

        delta_in_W = delta_in(graph, vertices)
        rpsw = infinity
        sigma: List = []

        if len(vertex_list) == 1 and delta_in_W <= k:
            rpsw = delta_in_W
            sigma = [vertex_list[0]]
        else:
            components = list(nx.weakly_connected_components(subgraph))

            if len(components) > 1:
                rpsw_list: List[int] = []
                for U_i in components:
                    rpsw_i, sigma_i = self._restricted_partial_scanwidth(
                        graph=graph,
                        vertices=U_i,
                        k=k,
                        infinity=infinity,
                    )
                    rpsw_list.append(rpsw_i)
                    sigma = sigma + sigma_i
                rpsw = max(rpsw_list)
            elif len(vertex_list) > 1 and delta_in_W <= k:
                for rho in roots:
                    new_vertices = vertex_list.copy()
                    new_vertices.remove(rho)
                    rpsw1, sigma1 = self._restricted_partial_scanwidth(
                        graph=graph,
                        vertices=new_vertices,
                        k=k,
                        infinity=infinity,
                    )
                    rpsw_prime = max(rpsw1, delta_in_W)
                    if rpsw_prime < rpsw:
                        rpsw = rpsw_prime
                        sigma = sigma1 + [rho]

        return rpsw, sigma
