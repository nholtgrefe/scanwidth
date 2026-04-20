"""Core DAG data object for the scanwidth package."""

from __future__ import annotations

from typing import Optional

import networkx as nx


class DAG:
    """Core directed acyclic graph data object.

    The :class:`DAG` is a lightweight wrapper around a ``networkx.DiGraph``
    intended purely as a data carrier. It exposes no algorithm methods;
    algorithms are implemented in dedicated solver classes in
    ``scanwidth.edge_scanwidth``.

    Parameters
    ----------
    graph : Optional[nx.DiGraph], optional
        NetworkX directed graph instance. If None, an empty directed graph
        is created.

    Raises
    ------
    TypeError
        If ``graph`` is not a ``networkx.DiGraph``.
    ValueError
        If ``graph`` is directed but contains at least one cycle.
    """

    def __init__(
        self,
        graph: Optional[nx.DiGraph] = None,
    ) -> None:
        if graph is None:
            self._graph = nx.DiGraph()
        else:
            if not isinstance(graph, nx.DiGraph):
                raise TypeError(
                    "graph must be a networkx.DiGraph or None."
                )
            if not nx.is_directed_acyclic_graph(graph):
                raise ValueError("graph must be a directed acyclic graph (DAG).")
            self._graph = graph

    @property
    def graph(self) -> nx.DiGraph:
        """Return the underlying ``networkx.DiGraph``.

        Returns
        -------
        nx.DiGraph
            Internal directed graph.
        """
        return self._graph

