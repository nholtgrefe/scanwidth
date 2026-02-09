"""Main module for computing scanwidth of directed acyclic graphs."""

import itertools
import os
import random
import sys
import time
from typing import Dict, List, Optional, Set, Tuple, Union

import networkx as nx
import numpy as np
from networkx.algorithms.flow import shortest_augmenting_path


# Global table for memoization
rpsw_table: Dict[str, Tuple[int, List]] = {}


class DAG:
    """Class for directed acyclic graphs to compute scanwidth.
    
    This class provides various algorithms and heuristics for computing
    the scanwidth of a DAG, including exact algorithms and approximation
    methods.
    """
    
    def __init__(
        self, 
        graph: Optional[Union[nx.DiGraph, str]] = None
    ) -> None:
        """Initialize a DAG object.
        
        Parameters
        ----------
        graph : Optional[Union[nx.DiGraph, str]], optional
            Either a NetworkX DiGraph object, a path to a text file in
            edge-list format (each line contains an arc of the graph),
            or None for an empty graph. Default is None.
        
        Examples
        --------
        >>> import networkx as nx
        >>> g = nx.DiGraph([(1, 2), (2, 3)])
        >>> dag = DAG(g)
        >>> dag = DAG("path/to/graph.el")
        >>> dag = DAG()  # Empty graph
        """
        if graph is None:  # Empty graph
            self.graph = nx.DiGraph()
        elif isinstance(graph, str):  # Initialize with graph from textfile
            self.graph = nx.DiGraph()
            edge_list = []
            with open(graph, 'r') as infile:
                data = infile.readlines()
                for i in data:
                    edge = i.split()
                    edge_list.append(tuple(edge))
            self.graph.add_edges_from(edge_list)
        else:  # Initialize with given graph
            self.graph = graph
        
        self.infinity = self.graph.number_of_edges() + 1
        self._memory = True
        self._opt_start_time = 0
        self.func_durations: List[float] = []
    
    def save_file(self, file_name: str) -> None:
        """Save the graph in edge-list format.
        
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
            for edge in list(self.graph.edges()):
                f.write(f"{edge[0]} {edge[1]}\n")
    
    def optimal_scanwidth(
        self, 
        reduced: bool = True, 
        method: int = 5, 
        memory: bool = True
    ) -> Tuple[Optional[int], Optional['Extension']]:
        """Compute the scanwidth of the DAG and return the corresponding extension.
        
        Parameters
        ----------
        reduced : bool, optional
            If True, first reduce the graph. Default is True.
        method : int, optional
            Algorithm method to use:
            - 1: Exhaustive search on all possible extensions
            - 2: Recursive algorithm (3-partition)
            - 3: Dynamic programming without component splitting
            - 4: Dynamic programming with component splitting
            - 5: Dynamic XP-program with increasing k
            Default is 5.
        memory : bool, optional
            If True, use extra space to increase speed. Default is True.
        
        Returns
        -------
        Tuple[Optional[int], Optional[Extension]]
            A tuple containing the scanwidth value and the corresponding
            Extension object. Returns (None, None) if computation fails.
        """
        if reduced:
            best_sw, best_ext = self._reduce(
                "optimal", 
                reduced=False, 
                method=method, 
                memory=memory
            )
            if best_sw is None:
                return None, None
            return best_sw, best_ext
        
        self._memory = memory
        global rpsw_table

        if method == 1:
            extension_list = [
                ext[::-1] for ext in nx.all_topological_sorts(self.graph)
            ]
            SW_sigma_list = []
            
            for sigma in extension_list:
                ext = Extension(self.graph, sigma)
                SW_sigma = ext.scanwidth()
                SW_sigma_list.append(SW_sigma)

            best_sw = min(SW_sigma_list)
            best_sigma_index = SW_sigma_list.index(best_sw)
            best_ext = Extension(self.graph, extension_list[best_sigma_index])
            
            return best_sw, best_ext
    
        elif method in [2, 3, 4]:
            rpsw_table = {}
            
            if method == 4:
                sw, sigma = self._restricted_partial_scanwidth(
                    self.graph.nodes(), 
                    cs=True
                )
            elif method == 2:
                sw, sigma = self._partial_scanwidth(
                    set(), 
                    set(self.graph.nodes()), 
                    set()
                )
            elif method == 3:
                sw, sigma = self._restricted_partial_scanwidth(
                    self.graph.nodes(), 
                    cs=False
                )
                        
            ext = Extension(self.graph, sigma)
            
            return sw, ext
            
        elif method == 5:
            rpsw_table = {}
                        
            for i in range(1, self.infinity + 1):
                res = self.optimal_k_scanwidth(i)
                    
                if res is False:
                    # Delete infinity values from table, they have to be recalculated.
                    for key in list(rpsw_table):
                        if rpsw_table[key][0] == self.infinity:
                            del rpsw_table[key]                 
                else:                    
                    sw, ext = res
                    return sw, ext
        
        return None, None
        
    def optimal_k_scanwidth(
        self, 
        k: int, 
        clear_table: bool = False
    ) -> Union[Tuple[int, 'Extension'], bool]:
        """Solve the fixed parameter version of scanwidth with rpsw_XP algorithm.
        
        Parameters
        ----------
        k : int
            The fixed parameter k for scanwidth.
        clear_table : bool, optional
            If True, clear the memoization table before computation.
            Must be True if calling this function by itself. Default is False.
        
        Returns
        -------
        Union[Tuple[int, Extension], bool]
            A tuple containing the scanwidth value and Extension object if
            scanwidth <= k, otherwise False.
        """
        if clear_table:
            global rpsw_table
            rpsw_table = {}
            
        sw, sigma = self._restricted_partial_scanwidth(self.graph.nodes(), k)
        
        if sw == self.infinity:
            return False
        else:
            ext = Extension(self.graph, sigma)
            return sw, ext
        
    def _restricted_partial_scanwidth(
        self, 
        vertices: Union[Set, List], 
        k: Optional[int] = None, 
        cs: bool = True
    ) -> Tuple[int, List]:
        """Compute restricted partial scanwidth.
        
        If k is None, outputs the rpsw and corresponding extension for the
        set vertices. Otherwise, only does this if rpsw <= k, otherwise
        it outputs 'infinity'.
        
        Parameters
        ----------
        vertices : Union[Set, List]
            Set of vertices to compute scanwidth for.
        k : Optional[int], optional
            Fixed parameter. If None, no fixed parameter is used.
            Default is None.
        cs : bool, optional
            If True, uses component splitting. Default is True.
        
        Returns
        -------
        Tuple[int, List]
            A tuple containing the scanwidth value and the extension (sigma).
        """
        # For the case that k = None, so no fixed parameter
        if k is None:
            k = self.infinity
        # For the fixed parameter version, we must use CS
        if k != self.infinity:
            assert cs is True
        
        vertex_list = list(vertices)  # Takes on the role of W
        subgraph = self.graph.subgraph(vertex_list)  # Takes on the role of G[W]      
        roots = [v for v in subgraph.nodes() if subgraph.in_degree(v) == 0]

        # Check table
        key = ""
        if self._memory:
            key = repr(sorted(roots))
            if key in rpsw_table.keys():
                return rpsw_table[key]

        delta_in_W = self._delta_in(vertices)
        
        # Initialize
        rpsw = self.infinity
        sigma: List = []
                
        # If |W| = 1 and delta_in(W) <= k
        if len(vertex_list) == 1 and delta_in_W <= k:
            rpsw = delta_in_W
            sigma = [vertex_list[0]]
        
        else:
            components = [
                comp for comp in nx.weakly_connected_components(subgraph)
            ]
            
            # If G[W] is weakly disconnected and CS is True
            if cs and len(components) > 1:
                rpsw_list = []
                for U_i in components:
                    rpsw_i, sigma_i = self._restricted_partial_scanwidth(
                        U_i, k, cs
                    )
                    rpsw_list.append(rpsw_i)
                    sigma = sigma + sigma_i
                    
                rpsw = max(rpsw_list)
                        
            # If |W| > 1 and delta_in(W) <= k
            elif len(vertex_list) > 1 and delta_in_W <= k:
                # Minimize over the roots of the subgraph
                for rho in roots:
                    new_vertices = vertex_list.copy()
                    new_vertices.remove(rho)

                    rpsw1, sigma1 = self._restricted_partial_scanwidth(
                        new_vertices, k, cs
                    )
                    
                    if cs:
                        rpsw2 = delta_in_W
                    else:
                        component = self._find_component(components, rho)
                        rpsw2 = self._delta_in(component)
                        
                    rpsw_prime = max(rpsw1, rpsw2)
                    
                    if rpsw_prime < rpsw:
                        rpsw = rpsw_prime
                        sigma = sigma1 + [rho]

        # Save to table
        if self._memory:
            rpsw_table[key] = (rpsw, sigma)

        return rpsw, sigma
    
    def _partial_scanwidth(
        self, 
        L: Set, 
        W: Set, 
        R: Set
    ) -> Tuple[int, List]:
        """Return the psw and corresponding extension for the ordered 3-partition L, W, R.
        
        Parameters
        ----------
        L : Set
            Left partition.
        W : Set
            Middle partition.
        R : Set
            Right partition.
        
        Returns
        -------
        Tuple[int, List]
            A tuple containing the partial scanwidth and extension.
        """
        # Initialize
        psw = self.infinity
        sigma: List = []

        # If |W| = 1
        if len(W) == 1:
            L_U_W = W.union(L)
            (w,) = W  # unpacks the unique element in W
            subgraph = self.graph.subgraph(L_U_W)  # Takes on the role of G[L U W]
            components = [
                comp for comp in nx.weakly_connected_components(subgraph)
            ]
            component = self._find_component(components, w)
                        
            psw = self._delta_in(component)
        
            sigma = [w]
        
        # If |W| > 1
        elif len(W) > 1:
            size = len(W) // 2
            for W_prime in map(set, itertools.combinations(W, size)):
                    
                # Check if W' sqsubseteq W
                check = len([
                    (u, v) for u in W_prime 
                    for v in W.difference(W_prime) 
                    if (u, v) in self.graph.edges()
                ])
                if check == 0:
                    psw_1, sigma_1 = self._partial_scanwidth(
                        L, W_prime, R.union(W.difference(W_prime))
                    )
                    psw_2, sigma_2 = self._partial_scanwidth(
                        L.union(W_prime), W.difference(W_prime), R
                    )
                    psw_prime = max(psw_1, psw_2)
                    if psw_prime < psw:
                        psw = psw_prime
                        sigma = sigma_1 + sigma_2
        
        return psw, sigma
        
    def cut_splitting_heuristic(
        self, 
        reduced: bool = True
    ) -> Tuple[int, 'Extension']:
        """Heuristic to compute the scanwidth using recursive cut-splitting.
        
        Parameters
        ----------
        reduced : bool, optional
            If True, first reduce the graph. Default is True.
        
        Returns
        -------
        Tuple[int, Extension]
            A tuple containing the scanwidth value and Extension object.
        """
        if reduced:
            sw, ext = self._reduce("cut", reduced=False)
            return sw, ext
        
        sigma = self._recursive_cut_splitting()
        ext = Extension(self.graph, sigma)
        sw = ext.scanwidth()
        
        return sw, ext
    
    def _recursive_cut_splitting(self) -> List:
        """Recursively split the graph at a DAG-cut and merge the rest.
        
        Returns
        -------
        List
            The extension (sigma) as a list of vertices.
        """
        edge_cut, source_set, sink_set = self._minimum_DAG_cut(non_trivial=True)

        # If there is no non-trivial cut anymore, just pick any ordering
        if edge_cut is None:
            return list(reversed(list(nx.topological_sort(self.graph))))
        
        subgraph_top = self.graph.subgraph(source_set)
        subgraph_down = self.graph.subgraph(sink_set)
        subgraph_top = subgraph_top.copy()
        subgraph_down = subgraph_down.copy()        
        
        # Add super leaf and super root
        super_leaf = "super_leaf_" + str(random.getrandbits(64))
        super_root = "super_root_" + str(random.getrandbits(64))
        subgraph_top.add_node(super_leaf, merged=True)
        subgraph_down.add_node(super_root, merged=True)
                
        # Add edges to super nodes (or add weight)
        for (u, v) in edge_cut:
            
            if 'weight' not in self.graph[u][v].keys():
                weight = 1
            else:
                weight = self.graph[u][v]['weight']
                        
            if subgraph_top.has_edge(u, super_leaf):
                subgraph_top[u][super_leaf]['weight'] += weight
            else:
                subgraph_top.add_edge(u, super_leaf, weight=weight)
                
            if subgraph_down.has_edge(super_root, v):
                subgraph_down[super_root][v]['weight'] += weight
            else:
                subgraph_down.add_edge(super_root, v, weight=weight)            
        
        top = DAG(subgraph_top)
        down = DAG(subgraph_down)

        # Recurse
        sigma_top = top._recursive_cut_splitting()
        sigma_down = down._recursive_cut_splitting()
        sigma_top.remove(super_leaf)
        sigma_down.remove(super_root)
                    
        return sigma_down + sigma_top
        
    def _minimum_DAG_cut(
        self, 
        non_trivial: bool = True
    ) -> Tuple[Optional[List], Optional[Set], Optional[Set]]:
        """Find the smallest DAG cut.
        
        Parameters
        ----------
        non_trivial : bool, optional
            If True, only looks for cuts such that both sides of the cut
            do not consist of a single merged node. Default is True.
        
        Returns
        -------
        Tuple[Optional[List], Optional[Set], Optional[Set]]
            A tuple containing the edges in the cut, and the node sets
            of both sides. Returns (None, None, None) if no cut is found.
        """
        aux_graph = self.graph.copy()
        
        # Give all edges weight of 1
        for (u, v) in aux_graph.edges():
            if 'weight' not in aux_graph[u][v].keys():
                aux_graph[u][v]['weight'] = 1
                
        # Add reverse edges with weight infinity
        aux_graph.add_weighted_edges_from(
            list((v, u, self.infinity) for (u, v) in aux_graph.edges)
        )

        roots = [v for v in self.graph.nodes() if self.graph.in_degree(v) == 0]
        leafs = [v for v in self.graph.nodes() if self.graph.out_degree(v) == 0]
        
        cut = self.infinity
        source_set: Set = set()
        sink_set: Set = set()
        size_diff = self.infinity + 1  # We prefer balanced separators if there is a tie
        
        # Find all DAG cuts
        if non_trivial is False:
            for root in roots:
                for leaf in leafs:
                    a, (b, c) = nx.minimum_cut(
                        aux_graph, 
                        root, 
                        leaf, 
                        capacity='weight',
                        flow_func=shortest_augmenting_path
                    )
                    d = abs(len(b) - len(c))
                    if a < cut or (a == cut and d < size_diff):
                        cut, source_set, sink_set, size_diff = a, b, c, d
                
        elif non_trivial is True:
            
            for root in roots:
                # If the root was merged, it is a trivial cut
                if 'merged' in self.graph.nodes[root].keys():
                    root_childs = list(self.graph.successors(root))
                else:
                    root_childs = [root]
                

                for leaf in leafs:
                    # If the leaf was merged, it is a trivial cut
                    if 'merged' in self.graph.nodes[leaf].keys():
                        leaf_parents = list(self.graph.predecessors(leaf))
                    else:
                        leaf_parents = [leaf]
                    
        
                    for root_child in root_childs:
                        for leaf_parent in leaf_parents:
                            
                            if root_child == leaf_parent:
                                continue  # Not a possible cut
                            a, (b, c) = nx.minimum_cut(
                                aux_graph, 
                                root_child, 
                                leaf_parent, 
                                capacity='weight', 
                                flow_func=shortest_augmenting_path
                            )
                            d = abs(len(b) - len(c))

                            if a < cut or (a == cut and d < size_diff):
                                cut, source_set, sink_set, size_diff = a, b, c, d
        
        if cut == self.infinity:
            # No DAG cut has been found
            return None, None, None
        
        # Find corresponding edges
        edge_cut = []
        for (u, v) in self.graph.edges():
            if u in source_set and v in sink_set:
                edge_cut.append((u, v))
        
        return edge_cut, source_set, sink_set
    
    def simulated_annealing(
        self, 
        max_iter: int = 100, 
        p_in: float = 0.9, 
        p_stop: float = 0.01, 
        init_ext: Union[str, 'Extension'] = 'greedy', 
        reduced: bool = False, 
        verbose: bool = True, 
        seed: int = 42
    ) -> Tuple[int, 'Extension', List[int]]:
        """Use simulated annealing heuristic to compute the scanwidth.
        
        Parameters
        ----------
        max_iter : int, optional
            Maximum number of iterations. Default is 100.
        p_in : float, optional
            The initial acceptance probability. Default is 0.9.
        p_stop : float, optional
            The final acceptance probability. Default is 0.01.
        init_ext : Union[str, Extension], optional
            Initial extension to be used. Can either be an Extension object,
            or a string ('cut', 'greedy' or 'random') to indicate which
            method should be used to get the initial extension. Default is 'greedy'.
        reduced : bool, optional
            Whether the graph should first be reduced. Default is False.
        verbose : bool, optional
            If True, print progress information. Default is True.
        seed : int, optional
            Random seed for reproducibility. Default is 42.
        
        Returns
        -------
        Tuple[int, Extension, List[int]]
            A tuple containing the best scanwidth value, the best Extension
            object, and a list of scanwidth values during the process.
        """
        if reduced:  # Split in blocks and do algorithm per block
            best_sw, best_ext = self._reduce(
                "annealing", 
                max_iter=max_iter, 
                p_in=p_in, 
                p_stop=p_stop, 
                init_ext=init_ext, 
                reduced=False, 
                verbose=verbose, 
                seed=seed
            )            
            return best_sw, best_ext, [best_sw]
        
        rng = np.random.RandomState(seed)
        
        if isinstance(init_ext, Extension):
            extension = init_ext
        elif init_ext == 'greedy':
            extension = self.greedy_heuristic()[1]
        elif init_ext == 'cut':
            extension = self.cut_splitting_heuristic()[1]
        elif init_ext == 'random':
            extension = self.random_extension()[1]
        else:
            raise ValueError(
                f"init_ext must be an Extension object or one of "
                f"{{'greedy', 'cut', 'random'}}, got {init_ext}"
            )
        
        tree_extension = extension.canonical_tree_extension()
        tree = tree_extension.tree
        
        # Initial scanwidth values
        sw_values: Dict = {}
        for v in self.graph.nodes:
            sw_values[v] = self._delta_in(nx.descendants(tree, v) | {v})
        
        # Algorithm initialization
        iteration_counter = 0
        best_tree = tree.copy()
        best_sw_values = sw_values.copy()
        best_sw = max(best_sw_values.values())
        epoch_length = len(self.graph.nodes())
        
        # Find out initial and stop temperature
        random_scanwidths = []
        for i in range(100):
            t = TreeExtension(self.graph, tree)
            s = t._random_neighbour_scanwidth(sw_values, rng)
            if s is not None: 
                random_scanwidths.append(max(s[0].values()))
                
        random_scanwidths = [
            elt - best_sw for elt in random_scanwidths if elt - best_sw > 0
        ]
        if len(random_scanwidths) == 0:
            delta_f = 1
        else:
            delta_f = sum(random_scanwidths) / len(random_scanwidths)
        init_temp = -delta_f / np.log(p_in)
        stop_temp = -delta_f / np.log(p_stop)
        temp = init_temp
        vals = [max(sw_values.values())]
        
        while temp >= stop_temp:
            
            # Per epoch, we do |V| iterations of the following
            for i in range(epoch_length):
                
                # Find a random neighbour and its scanwidth
                t = TreeExtension(self.graph, tree)
                res = t._random_neighbour_scanwidth(sw_values, rng)
                if res is None:
                    if verbose:
                        print("Only one possible tree extension")
                    return (
                        max(sw_values.values()), 
                        extension, 
                        [max(sw_values.values()) for i in range(max_iter + 1)]
                    )
                
                new_sw_values, vertex, parent, connected_vertices = res
                    
                # Calculate metropolis acceptance criterion
                diff = max(new_sw_values.values()) - max(sw_values.values())
                metropolis = np.exp(-diff / temp)
                
                # If improvement or random acceptance, accept the new tree
                if diff < 0 or rng.random() < metropolis:
                    
                    # Update values
                    sw_values = new_sw_values
                    
                    # Create new tree
                    succ = list(tree.successors(vertex))
                    grandparents = list(tree.predecessors(parent))
                    
                    for child in succ:
                        if child in connected_vertices:
                            tree.remove_edge(vertex, child)
                            tree.add_edge(parent, child)
                            
                    tree.remove_edge(parent, vertex)
                    
                    if len(grandparents) != 0:
                        grandparent = list(tree.predecessors(parent))[0]
                        tree.remove_edge(grandparent, parent)
                        tree.add_edge(grandparent, vertex)
                    
                    tree.add_edge(vertex, parent)
                    
                    # Keep track of best tree so far
                    if max(sw_values.values()) <= best_sw:
                        best_tree = tree.copy()
                        best_sw_values = sw_values.copy()
                        best_sw = max(best_sw_values.values())

                iteration_counter = iteration_counter + 1
                
                if verbose and iteration_counter % 100 == 0:
                    sys.stdout.write(
                        "\r" + f"current scanwidth: {max(sw_values.values())}; "
                        f"iteration {iteration_counter}; "
                        f"best scanwidth: {best_sw}"
                    )
                    time.sleep(1e-40)
            
            vals.append(max(sw_values.values()))
            
            # Change temperature
            alpha = (stop_temp / init_temp) ** (1 / (max_iter - 1))
            temp = alpha * temp            
        
        best_tree_extension = TreeExtension(self.graph, best_tree)

        best_extension = best_tree_extension.to_extension()
        return best_sw, best_extension, vals

    def greedy_heuristic(
        self, 
        reduced: bool = True
    ) -> Tuple[int, 'Extension']:
        """Greedy heuristic to find the scanwidth.
        
        Also returns the corresponding extension.
        
        Parameters
        ----------
        reduced : bool, optional
            If True, first reduce the graph. Default is True.
        
        Returns
        -------
        Tuple[int, Extension]
            A tuple containing the scanwidth value and Extension object.
        """
        if reduced:
            sw, ext = self._reduce("greedy", reduced=False)
            return sw, ext
        
        
        S = set(self.graph.nodes())
        T: Set = set()
        sigma: List = []
        sw = 0
        components_T: List[Set] = []  # Speedup: keeps track of the components of G[T]
        
        while len(S) > 0:
            x = None
            sw_x = self.infinity
            sub = self.graph.subgraph(S)
            leafs = [v for v in S if sub.out_degree(v) == 0]
            
            connection_dict: Dict = {}  # Keeps track of the components of G[T u {l}]
            for l in leafs:
                
                c = set(self.graph.successors(l))
                connected_vertices: Set = set()
                for comp in components_T:
                    if not c.isdisjoint(comp):  # If components has a child of l, it is connected
                        connected_vertices = connected_vertices | comp
                connected_vertices.add(l)
                connection_dict[l] = connected_vertices
                
                s = self._delta_in(connected_vertices)
                
                if s < sw_x:
                    x = l
                    sw_x = s
            
            S.remove(x)
            T.add(x)
            sigma = sigma + [x]
            sw = max(sw, sw_x)
            
            # New components of G[T]
            components_T = [
                comp for comp in components_T 
                if not comp.issubset(connection_dict[x])
            ]
            components_T.append(connection_dict[x])
        
        ext = Extension(self.graph, sigma)

        return sw, ext

    def random_extension(
        self, 
        seed: int = 42
    ) -> Tuple[int, 'Extension']:
        """Create a random extension of the DAG using the greedy structure.
        
        Also returns the corresponding scanwidth.
        
        Parameters
        ----------
        seed : int, optional
            Random seed for reproducibility. Default is 42.
        
        Returns
        -------
        Tuple[int, Extension]
            A tuple containing the scanwidth value and Extension object.
        """
        rng = np.random.RandomState(seed)
        
        S = set(self.graph.nodes())
        T: Set = set()
        sigma: List = []
        sw = 0
        components_T: List[Set] = []  # Speedup: keeps track of the components of G[T]
        
        while len(S) > 0:
            sub = self.graph.subgraph(S)
            leafs = [v for v in S if sub.out_degree(v) == 0]
            
            l = rng.choice(leafs)
                
            c = set(self.graph.successors(l))
            connected_vertices: Set = set()
            for comp in components_T:
                if not c.isdisjoint(comp):  # If components has a child of l, it is connected
                    connected_vertices = connected_vertices | comp
            connected_vertices.add(l)
            
            s = self._delta_in(connected_vertices)
            
            S.remove(l)
            T.add(l)
            sigma = sigma + [l]
            sw = max(sw, s)
            
            # New components of G[T]
            components_T = [
                comp for comp in components_T 
                if not comp.issubset(connected_vertices)
            ]
            components_T.append(connected_vertices)
        
        ext = Extension(self.graph, sigma)

        return sw, ext
    
    def sblocks(self) -> List[Set]:
        """Return a list of nodesets, one for each s-block of the graph.
        
        The list is ordered according to a reversed DFS in the sblock-cut-tree,
        starting from the rootblock. (Thus the rootblock is at the end.)
        This allows splitting the scanwidth problem over these blocks and then
        attaching them again in this ordering.
        
        Returns
        -------
        List[Set]
            A list of sets, where each set contains the nodes of an s-block.
        """
        roots = {v for v in self.graph.nodes() if self.graph.in_degree(v) == 0}
        aux = self.graph.to_undirected()

        # Create auxiliary graph, note that if G is rooted, this does not do anything.
        for root1 in roots:
            for root2 in roots:
                # Make a root-clique
                if (not (root1, root2) in aux.edges() and 
                    not (root2, root1) in aux.edges and 
                    root1 != root2):
                    aux.add_edge(root1, root2)
        
        sblock_sets = list(nx.biconnected_components(aux))
        dcut_vertices = list(nx.articulation_points(aux))
                        
        # Create sblock cut tree
        sblock_cut_tree = nx.Graph()
        for v in dcut_vertices:  # add the directed cut-vertices
            sblock_cut_tree.add_node(v)
        
        rootblock_index = None
        for i, block in enumerate(sblock_sets):  # add the s-blocks
            node_name = "block_" + str(i)
            if roots.issubset(block):  # keep track of the rootblock
                rootblock_index = i
            sblock_cut_tree.add_node(node_name)
            for v in dcut_vertices:
                if v in block:
                    sblock_cut_tree.add_edge(v, node_name)
        
        # DFS in the sblock cut tree
        sblock_order = list(
            nx.dfs_preorder_nodes(
                sblock_cut_tree, 
                source="block_" + str(rootblock_index)
            )
        )
        
        # Delete cut-vertices from order and only save indices:
        sblock_order = [
            int(name[6:]) for name in sblock_order 
            if str(name).startswith('block_')
        ]
        
        return [sblock_sets[i] for i in sblock_order]

    def _reduce(
        self, 
        func_type: str, 
        **kwargs
    ) -> Tuple[Optional[int], Optional['Extension']]:
        """Reduce the instance of the DAG, then apply scanwidth algorithm.
        
        This function first reduces the instance of the DAG, then it applies
        any scanwidth-algorithm or -heuristic on the reduced instance(s).
        
        Parameters
        ----------
        func_type : str
            Type of function to apply: 'optimal', 'cut', 'greedy', or 'annealing'.
        **kwargs
            Additional keyword arguments to pass to the function.
        
        Returns
        -------
        Tuple[Optional[int], Optional[Extension]]
            A tuple containing the scanwidth value and Extension object.
            Returns (None, None) if computation fails.
        """
        sblock_sets = self.sblocks()
        sigma: List = []
        sw = 0
        
        # We find the scanwidth on each block separately
        for sblock_set in sblock_sets:
            # This is a single edge, so we already know sw = 1
            if len(sblock_set) == 2: 
                u, v = sblock_set
                if (u, v) in self.graph.edges():
                    partial_sigma = [v, u]
                else:
                    partial_sigma = [u, v]
                sw = max(sw, 1)
            
            else:
                subgraph = self.graph.subgraph(sblock_set)
                roots = [
                    v for v in subgraph.nodes() 
                    if subgraph.in_degree(v) == 0 and subgraph.out_degree(v) == 2
                ]
                leafs = [
                    v for v in subgraph.nodes() 
                    if subgraph.out_degree(v) == 0 and subgraph.in_degree(v) == 2
                ]
                flow_nodes = [
                    v for v in subgraph.nodes() 
                    if subgraph.out_degree(v) == 1 and subgraph.in_degree(v) == 1
                ]
                
                # This means that the block is a directed cycle, so it has sw = 2
                if (len(roots) == 1 and len(leafs) == 1 and 
                    len(flow_nodes) == len(subgraph.nodes()) - 2):
                    partial_sigma = list(
                        reversed(list(nx.topological_sort(subgraph)))
                    )  # any ordering is okay
                    sw = max(sw, 2) 
                
                # This means the block is neither an edge or a cycle, so we apply the algorithm
                else: 
                    # First contract flow-edges and save the contractions
                    history = []
                    subgraph = subgraph.copy()
                    for v in flow_nodes:
                        u = list(subgraph.predecessors(v))[0]
                        w = list(subgraph.successors(v))[0]
                        if (u, w) not in subgraph.edges():
                            subgraph.remove_node(v)
                            subgraph.add_edge(u, w)
                            history.append((w, v))
                    
                    # Run algorithm to find sigma and sw                    
                    S = DAG(subgraph)

                    if func_type == 'optimal':
                        res = S.optimal_scanwidth(**kwargs)
                    elif func_type == 'cut':
                        res = S.cut_splitting_heuristic(**kwargs)
                    elif func_type == 'greedy':
                        res = S.greedy_heuristic(**kwargs)
                    elif func_type == 'annealing':
                        res = S.simulated_annealing(**kwargs)
                        if res is not None and len(res) == 3:
                            res = (res[0], res[1])
                    else:
                        raise ValueError(
                            f"Unknown func_type: {func_type}. "
                            f"Must be one of {{'optimal', 'cut', 'greedy', 'annealing'}}"
                        )
                    
                    if res is None or res[0] is None:
                        return None, None
                    block_sw, partial_ext = res
                    partial_sigma = partial_ext.sigma
                    
                    sw = max(sw, block_sw)
                    
                    # Backtrack the edge-contractions to get the partial extension for the whole block
                    history.reverse()
                    for (w, v) in history:
                        i = partial_sigma.index(w)
                        partial_sigma.insert(i + 1, v)
            
            # Add partial extension to the whole extension
            partial_sigma = [v for v in partial_sigma if v not in sigma]
            sigma = partial_sigma + sigma
            
        return sw, Extension(self.graph, sigma)

    def nr_of_leaves(self) -> int:
        """Return the number of leaves of the graph.
        
        Returns
        -------
        int
            The number of leaves (nodes with out-degree 0).
        """
        return len([u for u in self.graph.nodes() if self.graph.out_degree(u) == 0])

    def _delta_in(
        self, 
        vertex_set: Union[Set, List], 
        sink: bool = True
    ) -> int:
        """Return the indegree of vertex_set.
        
        Setting sink to True, if we already know that it is a sinkset
        will speed up computation.
        
        Parameters
        ----------
        vertex_set : Union[Set, List]
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

    @staticmethod
    def _find_component(components: List[Set], v) -> Set:
        """Return the component in the list components that contains v.
        
        Parameters
        ----------
        components : List[Set]
            List of component sets.
        v
            Vertex to find.
        
        Returns
        -------
        Set
            The component set containing v.
        """
        for comp in components:
            if v in comp:
                return comp
        return set()


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
