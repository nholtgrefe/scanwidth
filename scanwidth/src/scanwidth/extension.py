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
        """
        self.graph = graph
        
        if isinstance(sigma, str):  # Initialize with sigma from textfile
            self.sigma: List = []
            with open(sigma, 'r') as infile:
                data = infile.readlines()
                for i in data:
                    v = i.split()[0]
                    self.sigma.append(v)
                            
        else:  # Initialize with given sigma
            self.sigma = sigma
            
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
            for v in self.sigma:
                f.write(f"{v}\n")
        
    def scanwidth(self) -> int:
        """Calculate the scanwidth of the extension sigma for the DAG.
        
        Returns
        -------
        int
            The scanwidth value.
        """
        SW_i_list = []
        
        for i in range(len(self.sigma)):
            SW_i = self.scanwidth_at_vertex_i(i)
            SW_i_list.append(SW_i)
        
        return max(SW_i_list)

    def scanwidth_at_vertex_i(
        self, 
        i: int, 
        position: bool = True
    ) -> int:
        """Calculate the size of the set SW of the extension sigma at position i.
        
        If position is False, we find the scanwidth at the vertex i
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
            The scanwidth at position/vertex i.
        """
        if not position:
            i = self.sigma.index(i)
        
        left = self.sigma[0:i + 1]
        
        sub = self.graph.subgraph(left)
        components = [comp for comp in nx.weakly_connected_components(sub)]
        connected_vertices = set()
        for comp in components:
            if self.sigma[i] in comp:
                connected_vertices = comp
                break
                
        SW_i = 0
        for w in connected_vertices:
            SW_i = SW_i + self.graph.in_degree(w) - self.graph.out_degree(w)
        
        return SW_i
    
    def canonical_tree_extension(self) -> TreeExtension:
        """Create the canonical tree extension with the same scanwidth as sigma.
        
        Returns
        -------
        TreeExtension
            A TreeExtension object with the same scanwidth as this extension.
        """
        # Initialize
        Gamma = nx.DiGraph()
        sig = self.sigma.copy()
        rho: Dict = {node: None for node in self.graph.nodes()}
        
        while len(sig) > 0:
            v = sig[0]
            sig.remove(v)
            C = list(self.graph.successors(v))
            Gamma.add_node(v)
            rho[v] = v
            
            if len(C) > 0:
                R = set([rho[c] for c in C])
                for r in R:
                    Gamma.add_edge(v, r)
                for u in Gamma.nodes():
                    if rho[u] in R:
                        rho[u] = v
        
        tree = TreeExtension(self.graph, Gamma)
        
        return tree
    
    def is_extension(self) -> bool:
        """Check if sigma is indeed an extension of the graph.
        
        Returns
        -------
        bool
            True if sigma is a valid extension, False otherwise.
        """
        seen = []
        for i in range(len(self.sigma)):
            v = self.sigma[i]
            succ = list(self.graph.successors(v))
            for w in succ:
                if w not in seen:
                    return False
            seen.append(v)
        return True
