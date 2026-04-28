"""Extension class for extensions of graphs."""

from typing import Dict, List

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.tree_extension import TreeExtension


class Extension:
    """Class for an extension of a graph, starting with the leafnodes."""
    
    def __init__(
        self, 
        dag: DAG, 
        ordering: List
    ) -> None:
        """Initialize an Extension object.
        
        Parameters
        ----------
        dag : DAG
            The corresponding DAG wrapper object.
        ordering : List
            List of nodes representing the extension order.

        Raises
        ------
        TypeError
            If ``dag`` is not a :class:`DAG`.
        TypeError
            If ``ordering`` is not a list.
        ValueError
            If ``ordering`` is not a valid extension of ``dag``.
        """
        if not isinstance(dag, DAG):
            raise TypeError("dag must be a DAG instance.")
        if not isinstance(ordering, list):
            raise TypeError("ordering must be a list.")
        self._dag = dag
        self._ordering = ordering

        if not self._is_valid():
            raise ValueError("ordering is not a valid extension of graph.")

    @property
    def dag(self) -> DAG:
        """Return the underlying DAG object."""
        return self._dag

    @property
    def ordering(self) -> List:
        """Return the extension order."""
        return self._ordering.copy()
            
    def edge_scanwidth(self) -> int:
        """Calculate edge scanwidth of the extension ordering for the DAG.
        
        Returns
        -------
        int
            The scanwidth value.
        """
        return max(len(self.edge_scanwidth_bag(v)) for v in self._ordering)

    def node_scanwidth(self) -> int:
        """Calculate node scanwidth of the extension ordering for the DAG.

        Returns
        -------
        int
            The node scanwidth value.
        """
        return max(len(self.node_scanwidth_bag(v)) for v in self._ordering)

    def edge_scanwidth_bag(self, vertex: object) -> set:
        """Return the set of edges in the edge scanwidth bag for a vertex.

        Parameters
        ----------
        vertex : object
            Vertex in the extension order ``ordering``.

        Returns
        -------
        set
            Set of edges in ``SW_v``.

        Raises
        ------
        ValueError
            If ``vertex`` is not in the extension order.
        """
        if vertex not in self._ordering:
            raise ValueError("vertex must be a node in the extension order.")

        connected_vertices = self._connected_vertices_for(vertex)

        return {
            (u, w)
            for (u, w) in self._dag.graph.edges()
            if u not in connected_vertices and w in connected_vertices
        }

    def node_scanwidth_bag(self, vertex: object) -> set:
        """Return node scanwidth bag for a vertex in the extension ordering.

        Parameters
        ----------
        vertex : object
            Vertex in the extension order ``ordering``.

        Returns
        -------
        set
            Set of parent vertices of edges in ``SW_v``.

        Raises
        ------
        ValueError
            If ``vertex`` is not in the extension order.
        """
        if vertex not in self._ordering:
            raise ValueError("vertex must be a node in the extension order.")

        connected_vertices = self._connected_vertices_for(vertex)

        return {
            u
            for (u, w) in self._dag.graph.edges()
            if u not in connected_vertices and w in connected_vertices
        }

    def _connected_vertices_for(self, vertex: object) -> set:
        """Return connected component of ``vertex`` in current ordering prefix."""
        i = self._ordering.index(vertex)
        left = self._ordering[0:i + 1]

        sub = self._dag.graph.subgraph(left)
        components = [comp for comp in nx.weakly_connected_components(sub)]
        for comp in components:
            if vertex in comp:
                return comp
        return set()
    
    def to_canonical_tree_extension(self) -> TreeExtension:
        """Create canonical tree extension with same edge scanwidth as ordering.
        
        Returns
        -------
        TreeExtension
            A TreeExtension object with the same scanwidth as this extension.
        """
        # Initialize
        Gamma = nx.DiGraph()
        sig = self._ordering.copy()
        rho: Dict = {node: None for node in self._dag.graph.nodes()}
        
        while len(sig) > 0:
            v = sig[0]
            sig.remove(v)
            C = list(self._dag.graph.successors(v))
            Gamma.add_node(v)
            rho[v] = v
            
            if len(C) > 0:
                R = set([rho[c] for c in C])
                for r in R:
                    Gamma.add_edge(v, r)
                for u in Gamma.nodes():
                    if rho[u] in R:
                        rho[u] = v
        
        tree = TreeExtension(self._dag, Gamma)
        
        return tree

    def canonical_tree_extension(self) -> TreeExtension:
        """Return the canonical tree extension.

        This forwards to :meth:`to_canonical_tree_extension`.
        """
        return self.to_canonical_tree_extension()
    
    def _is_valid(self) -> bool:
        """Check if ordering is indeed an extension of the graph.
        
        Returns
        -------
        bool
            True if ordering is a valid extension, False otherwise.
        """
        seen = []
        if len(self._ordering) != len(self._dag.graph.nodes()):
            return False
        if set(self._ordering) != set(self._dag.graph.nodes()):
            return False

        for i in range(len(self._ordering)):
            v = self._ordering[i]
            succ = list(self._dag.graph.successors(v))
            for w in succ:
                if w not in seen:
                    return False
            seen.append(v)
        return True

