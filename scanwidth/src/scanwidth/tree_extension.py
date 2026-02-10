"""TreeExtension class for tree extensions of graphs."""

from typing import Dict, Optional, Set, Tuple, TYPE_CHECKING

import networkx as nx
import numpy as np

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
        """
        self.graph = graph
        self.tree = tree
    
    def scanwidth(self) -> int:
        """Calculate the scanwidth of the tree extension for the DAG.
        
        Returns
        -------
        int
            The scanwidth value.
        """
        GW_v_list = []
        
        for v in self.tree.nodes():
            GW_v = self.scanwidth_at_vertex(v)
            GW_v_list.append(GW_v)
        
        return max(GW_v_list)

    def scanwidth_at_vertex(self, v) -> int:
        """Calculate the scanwidth of the tree extension at vertex v.
        
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
    
    def _random_neighbour_scanwidth(
        self, 
        sw_values: Dict, 
        rng: np.random.RandomState
    ) -> Optional[Tuple[Dict, object, object, Set]]:
        """Find the scanwidth of a random neighbour of the tree extension.
        
        Takes as input a dictionary of the scanwidth values of the current tree
        and a random number generator rng. Returns a dictionary of new
        scanwidth-values, the two vertices that are swapped and a set of vertices
        that are connected to the parent in the graph G[1..parent].
        
        Parameters
        ----------
        sw_values : Dict
            Dictionary of scanwidth values for the current tree.
        rng : np.random.RandomState
            Random number generator.
        
        Returns
        -------
        Optional[Tuple[Dict, object, object, Set]]
            A tuple containing new scanwidth values dictionary, vertex, parent,
            and connected_vertices set. Returns None if no valid neighbour exists.
        """
        # Choose random vertex
        possible_choices = []
        for v in self.graph.nodes:
            pred = list(self.tree.predecessors(v))
            if len(pred) == 0:  # v is the root and we can not move the root up
                continue
            elif (pred[0], v) in self.graph.edges:  # If an edge from pred to vertex, we can not swap them
                continue
            else:
                possible_choices.append(v)
                
        if len(possible_choices) == 0:
            return None
        
        # Vertex to be swapped above
        vertex = rng.choice(possible_choices)
        parent = list(self.tree.predecessors(vertex))[0]
        
        # Calculate the new sw_values (only parent and vertex change)
        new_sw_values = sw_values.copy()
        
        # New value for vertex
        new_sw_values[vertex] = sw_values[parent]
        
        # New value for parent
        desc = list(nx.descendants(self.tree, parent))
        desc.remove(vertex)
        desc.append(parent)
        
        connected_vertices = {parent}
        for child in self.tree.successors(parent):
            if child != vertex:
                connected_vertices = (
                    connected_vertices | nx.descendants(self.tree, child) | {child}
                )
        for child in self.tree.successors(vertex):
            sinkset = nx.descendants(self.tree, child) | {child}
            connected = False
            
            for u in sinkset:
                if (parent, u) in self.graph.edges():
                    connected = True
                    break
            
            if connected:
                connected_vertices = (
                    connected_vertices | nx.descendants(self.tree, child) | {child}
                )

        new_sw_values[parent] = self._delta_in(connected_vertices)
        
        return new_sw_values, vertex, parent, connected_vertices
        
    def _delta_in(
        self, 
        vertex_set: Set, 
        sink: bool = True
    ) -> int:
        """Return the indegree of vertex_set.
        
        Setting sink to True, if we already know that it is a sinkset
        will speed up computation.
        
        Parameters
        ----------
        vertex_set : Set
            Set of vertices to compute indegree for.
        sink : bool, optional
            If True, assumes vertex_set is a sinkset for optimization.
            Default is True.
        
        Returns
        -------
        int
            The indegree of the vertex set.
        """
        res = 0
        if sink:
            if len(vertex_set) < len(self.graph.nodes()) / 2:  # Return indegree of W
                for v in vertex_set:
                    res = res + self.graph.in_degree(v) - self.graph.out_degree(v)
            else:  # Return outdegree of V / W
                for v in self.graph.nodes:
                    if v not in vertex_set:
                        res = res - self.graph.in_degree(v) + self.graph.out_degree(v)
            
        if not sink:  # If no sinkset
            for (u, v) in self.graph.edges():
                if u not in vertex_set and v in vertex_set:
                    res = res + 1            
            
        return res
