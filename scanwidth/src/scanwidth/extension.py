"""Extension class for extensions of graphs."""

import os
from typing import Dict, List, Union

import networkx as nx

from scanwidth.tree_extension import TreeExtension


class Extension:
    """Class for an extension of a graph, starting with the leafnodes."""
    
    def __init__(
        self, 
        graph: nx.DiGraph, 
        sigma: Union[List, str]
    ) -> None:
        """Initialize an Extension object.
        
        Parameters
        ----------
        graph : nx.DiGraph
            The corresponding graph as a NetworkX DiGraph.
        sigma : Union[List, str]
            A list of nodes for the extension, or a path to a text file
            containing the extension (one vertex per line).

        Raises
        ------
        ValueError
            If ``sigma`` is not a valid extension of ``graph``.
        """
        self._graph = graph
        
        if isinstance(sigma, str):  # Initialize with sigma from textfile
            self._sigma: List = []
            with open(sigma, 'r') as infile:
                data = infile.readlines()
                for i in data:
                    v = i.split()[0]
                    self._sigma.append(v)
                            
        else:  # Initialize with given sigma
            self._sigma = sigma

        if not self._is_extension():
            raise ValueError("sigma is not a valid extension of graph.")

    @property
    def graph(self) -> nx.DiGraph:
        """Return the underlying graph."""
        return self._graph

    @property
    def sigma(self) -> List:
        """Return the extension order."""
        return self._sigma.copy()
            
    def save_file(self, file_name: str) -> None:
        """Save the extension in a file.
        
        Each line contains one vertex of the extension.
        
        Parameters
        ----------
        file_name : str
            Path where the file should be saved.
        
        Raises
        ------
        ValueError
            If the file already exists.
        """
        if os.path.exists(file_name):
            raise ValueError("File already exists.")
        
        with open(file_name, "w+") as f:
            for v in self._sigma:
                f.write(f"{v}\n")
        
    def edge_scanwidth(self) -> int:
        """Calculate edge scanwidth of the extension sigma for the DAG.
        
        Returns
        -------
        int
            The scanwidth value.
        """
        SW_i_list = []
        
        for i in range(len(self._sigma)):
            SW_i = self.edge_scanwidth_at_vertex_i(i)
            SW_i_list.append(SW_i)
        
        return max(SW_i_list)

    def edge_scanwidth_at_vertex_i(
        self, 
        i: int, 
        position: bool = True
    ) -> int:
        """Calculate edge scanwidth at position (or vertex) ``i`` in sigma.
        
        If position is False, we find edge scanwidth at the vertex i
        (i.e. the node-name).
        
        Parameters
        ----------
        i : int
            Position index or vertex name (if position=False).
        position : bool, optional
            If True, i is treated as a position index. If False, i is
            treated as a vertex name. Default is True.
        
        Returns
        -------
        int
            The edge scanwidth at position/vertex i.
        """
        if not position:
            i = self._sigma.index(i)
        
        left = self._sigma[0:i + 1]
        
        sub = self._graph.subgraph(left)
        components = [comp for comp in nx.weakly_connected_components(sub)]
        connected_vertices = set()
        for comp in components:
            if self._sigma[i] in comp:
                connected_vertices = comp
                break
                
        SW_i = 0
        for w in connected_vertices:
            SW_i = SW_i + self._graph.in_degree(w) - self._graph.out_degree(w)
        
        return SW_i
    
    def canonical_tree_extension(self) -> TreeExtension:
        """Create canonical tree extension with same edge scanwidth as sigma.
        
        Returns
        -------
        TreeExtension
            A TreeExtension object with the same scanwidth as this extension.
        """
        # Initialize
        Gamma = nx.DiGraph()
        sig = self._sigma.copy()
        rho: Dict = {node: None for node in self._graph.nodes()}
        
        while len(sig) > 0:
            v = sig[0]
            sig.remove(v)
            C = list(self._graph.successors(v))
            Gamma.add_node(v)
            rho[v] = v
            
            if len(C) > 0:
                R = set([rho[c] for c in C])
                for r in R:
                    Gamma.add_edge(v, r)
                for u in Gamma.nodes():
                    if rho[u] in R:
                        rho[u] = v
        
        tree = TreeExtension(self._graph, Gamma)
        
        return tree
    
    def _is_extension(self) -> bool:
        """Check if sigma is indeed an extension of the graph.
        
        Returns
        -------
        bool
            True if sigma is a valid extension, False otherwise.
        """
        seen = []
        if len(self._sigma) != len(self._graph.nodes()):
            return False
        if set(self._sigma) != set(self._graph.nodes()):
            return False

        for i in range(len(self._sigma)):
            v = self._sigma[i]
            succ = list(self._graph.successors(v))
            for w in succ:
                if w not in seen:
                    return False
            seen.append(v)
        return True

