"""TreeExtension class for tree extensions of graphs."""

from typing import TYPE_CHECKING

import networkx as nx

from scanwidth.dag import DAG

if TYPE_CHECKING:
    from scanwidth.extension import Extension


class TreeExtension:
    """Class for a tree extension of a graph."""
    
    def __init__(
        self, 
        dag: DAG, 
        tree: nx.DiGraph
    ) -> None:
        """Initialize a TreeExtension object.
        
        Parameters
        ----------
        dag : DAG
            The corresponding DAG wrapper object.
        tree : nx.DiGraph
            The tree extension as a NetworkX DiGraph.

        Raises
        ------
        TypeError
            If ``dag`` is not a :class:`DAG`.
        ValueError
            If ``tree`` is not a valid tree extension of ``dag``.
        """
        if not isinstance(dag, DAG):
            raise TypeError("dag must be a DAG instance.")
        self._dag = dag
        self._tree = tree
        if not self._is_valid():
            raise ValueError("tree is not a valid tree extension of graph.")

    @property
    def dag(self) -> DAG:
        """Return the underlying DAG object."""
        return self._dag

    @property
    def tree(self) -> nx.DiGraph:
        """Return the tree extension graph."""
        return self._tree
    
    def edge_scanwidth(self) -> int:
        """Calculate the edge scanwidth of the tree extension for the DAG.
        
        Returns
        -------
        int
            The scanwidth value.
        """

        return max(len(self.edge_scanwidth_bag(v)) for v in self.tree.nodes())

    def node_scanwidth(self) -> int:
        """Calculate node scanwidth of the tree extension for the DAG.

        Returns
        -------
        int
            The node scanwidth value.
        """
        return max(len(self.node_scanwidth_bag(v)) for v in self.tree.nodes())

    def edge_scanwidth_bag(self, vertex: object) -> set:
        """Return the set of edges in the edge scanwidth bag for a tree vertex.
        
        Parameters
        ----------
        vertex : object
            Vertex of the tree extension.
        
        Returns
        -------
        set
            Set of edges in ``GW_v``.

        Raises
        ------
        ValueError
            If ``vertex`` is not a vertex of ``self.tree``.
        """
        if vertex not in self.tree.nodes:
            raise ValueError("vertex must be a node of the tree extension.")

        left = nx.descendants(self.tree, vertex)
        left.add(vertex)
        return {
            (u, w)
            for (u, w) in self.dag.graph.edges()
            if u not in left and w in left
        }

    def node_scanwidth_bag(self, vertex: object) -> set:
        """Return node scanwidth bag for a tree extension vertex.

        Parameters
        ----------
        vertex : object
            Vertex of the tree extension.

        Returns
        -------
        set
            Set of parent vertices of edges in ``GW_v``.

        Raises
        ------
        ValueError
            If ``vertex`` is not a vertex of ``self.tree``.
        """
        if vertex not in self.tree.nodes:
            raise ValueError("vertex must be a node of the tree extension.")

        left = nx.descendants(self.tree, vertex)
        left.add(vertex)
        return {
            u
            for (u, w) in self.dag.graph.edges()
            if u not in left and w in left
        }
    
    def to_extension(self) -> 'Extension':
        """Return an extension of the tree.
        
        If the tree extension is canonical, this extension gives the same scanwidth.
        
        Returns
        -------
        Extension
            An Extension object corresponding to this tree extension.
        """
        from scanwidth.extension import Extension
        
        sigma = list(reversed(list(nx.topological_sort(self.tree))))
        return Extension(self.dag, sigma)
    
    def is_canonical(self) -> bool:
        """Report whether the tree extension is canonical.
        
        Returns
        -------
        bool
            True if the tree extension is canonical, False otherwise.
        """
        for v in self.tree.nodes():
            left = nx.descendants(self.tree, v)
            left.add(v)
            if not nx.is_weakly_connected(self.dag.graph.subgraph(left)):
                return False
            
        return True

    def _is_valid(self) -> bool:
        """Return whether ``self.tree`` is a valid tree extension of ``graph``."""
        graph_nodes = set(self.dag.graph.nodes())
        tree_nodes = set(self.tree.nodes())
        if graph_nodes != tree_nodes:
            return False

        if len(tree_nodes) == 0:
            return True

        if not nx.is_directed_acyclic_graph(self.tree):
            return False
        if self.tree.number_of_edges() != len(tree_nodes) - 1:
            return False
        if not nx.is_weakly_connected(self.tree):
            return False

        roots = [v for v in self.tree.nodes() if self.tree.in_degree(v) == 0]
        if len(roots) != 1:
            return False
        for v in self.tree.nodes():
            indeg = self.tree.in_degree(v)
            if v == roots[0]:
                if indeg != 0:
                    return False
            elif indeg != 1:
                return False

        for (u, v) in self.dag.graph.edges():
            if not nx.has_path(self.tree, u, v):
                return False
        return True
    
