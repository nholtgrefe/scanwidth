"""TreeExtension class for tree extensions of graphs."""

from typing import TYPE_CHECKING

import networkx as nx

if TYPE_CHECKING:
    from scanwidth.extension import Extension


class TreeExtension:
    """Class for a tree extension of a graph."""
    
    def __init__(
        self, 
        graph: nx.DiGraph, 
        tree: nx.DiGraph
    ) -> None:
        """Initialize a TreeExtension object.
        
        Parameters
        ----------
        graph : nx.DiGraph
            The corresponding graph as a NetworkX DiGraph.
        tree : nx.DiGraph
            The tree extension as a NetworkX DiGraph.

        Raises
        ------
        ValueError
            If ``tree`` is not a valid tree extension of ``graph``.
        """
        self._graph = graph
        self._tree = tree
        if not self._is_tree_extension():
            raise ValueError("tree is not a valid tree extension of graph.")

    @property
    def graph(self) -> nx.DiGraph:
        """Return the underlying graph."""
        return self._graph

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
        GW_v_list = []
        
        for v in self.tree.nodes():
            GW_v = self.edge_scanwidth_at_vertex(v)
            GW_v_list.append(GW_v)
        
        return max(GW_v_list)

    def edge_scanwidth_at_vertex(self, v) -> int:
        """Calculate edge scanwidth of the tree extension at vertex ``v``.
        
        Parameters
        ----------
        v
            The vertex to compute scanwidth at.
        
        Returns
        -------
        int
            The scanwidth at vertex v.
        """
        left = nx.descendants(self.tree, v)
        left.add(v)
        
        GW_v = 0
        for w in left:
            GW_v = GW_v + self.graph.in_degree(w) - self.graph.out_degree(w)
            
        return GW_v
    
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
        return Extension(self.graph, sigma)
    
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
            if not nx.is_weakly_connected(self.graph.subgraph(left)):
                return False
            
        return True

    def _is_tree_extension(self) -> bool:
        """Return whether ``self.tree`` is a valid tree extension of ``graph``."""
        graph_nodes = set(self.graph.nodes())
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

        for (u, v) in self.graph.edges():
            if not nx.has_path(self.tree, u, v):
                return False
        return True
    
