import random
import networkx as nx
import numpy as np
from scanwidth import Network

'''
Code for generating birth-hybridization network (ZODS-generator)
- Adapted from Esther Julien's code for the paper:
    van Iersel, L., Jones, M., Julien, E., & Murakami, Y. (2023). Making a Network Orchard by Adding Leaves. 
    arXiv preprint arXiv:2305.03106.
- She adapted the code from Remie Janssen's code for the paper (with help of Celine Scornavacca for distances on arcs):
    Janssen, R., & Liu, P. (2021). Comparing the topology of phylogenetic network generators. 
    Journal of Bioinformatics and Computational Biology, 19(06), 2140012.
- Originally, the method is from the paper:
    Zhang, C., Ogilvie, H.A., Drummond, A.J., Stadler, T.: Bayesian inference of species networks from multilocus
    sequence data. Molecular biology and evolution 35(2), 504–517 (2018)
'''


def ZODS_generator(time_limit, speciation_rate, hybridization_rate, rng, taxa_goal=None, max_retics=None):
    """ZODS-generator according to their birth-hybridization process. Parameters:
        - time_limit: Starting time of the process.
        - speciation_rate: how often species split into two.
        - hybridization_rate: how often two species merge.
        - rng: a random number generator.
        - taxa_goal: number of wanted taxa, i.e. leaves.
        - max_retics: maximum number of reticulations.
        Returns the network, number of reticulations and number of leaves."""
            
    nw = nx.DiGraph()
    nw.add_node(0)
    leaves = {0}
    current_node = 1
    no_of_leaves = 1
    retics = 0
    extra_time = rng.exponential(1/float(speciation_rate))
    current_time = extra_time
    current_speciation_rate = float(speciation_rate)
    current_hybridization_rate = hybridization_rate
    rate = current_speciation_rate + current_hybridization_rate

    while current_time < time_limit and (taxa_goal is not None and taxa_goal != no_of_leaves):
        # Speciate
        if (0 in leaves) or (rng.random() < current_speciation_rate / rate): 
            splitting_leaf = random.choice(list(leaves))
            nw.add_edges_from([(splitting_leaf, current_node), (splitting_leaf, current_node + 1)])
            leaves.remove(splitting_leaf)
            leaves.add(current_node)
            leaves.add(current_node+1)
            current_node += 2
            no_of_leaves += 1
            
        # Hybridize    
        elif len(leaves) >= 2: 
            merging = rng.choice(tuple(leaves), 2, replace=False)
            l0 = merging[0]
            l1 = merging[1]
            pl0 = -1
            for p in nw.predecessors(l0):
                pl0 = p
            pl1 = -1
            for p in nw.predecessors(l1):
                pl1 = p
            # If pl0==pl1, the new hybridization results in parallel edges.
            if pl0 != pl1:
                nw.remove_node(l1)
                nw.add_edges_from([(pl1, l0), (l0, current_node)])
                leaves.remove(l0)
                leaves.remove(l1)
                leaves.add(current_node)
                current_node += 1
                no_of_leaves -= 1
                retics += 1
        else:
            break
                
        current_speciation_rate = float(speciation_rate*no_of_leaves)
        current_hybridization_rate = float(hybridization_rate * (no_of_leaves * (no_of_leaves - 1))/2)
        rate = current_speciation_rate + current_hybridization_rate
        extra_time = rng.exponential(1/rate)
        current_time += extra_time
        if max_retics is not None and retics > max_retics:
            return None, retics, no_of_leaves

    # nothing has happened yet, and there is only one node
    if len(nw) == 1:
        nw.add_edges_from([(0, 1)])
        
    return Network(nw), retics, no_of_leaves


def create_random_network(nr_leaves, nr_reticulations, seed=42):
    """Iterates the ZODS-generator until it finds a network with the specified
    number of leaves and reticulations. Uses a speciation rate of 1, and a 
    hybridization rate in the interval [0.0001, 0.4]"""
    s_rate = 1.0 # We set speciation_rate to 1, else we just get indegree-outdegree-1 edges.
    rng = np.random.RandomState(seed) # Only sets seed for this generator
    
    while True:
        h_rate = rng.uniform(0.0001, 0.4)
        nw, rets, leaves = ZODS_generator(500, s_rate, h_rate, rng, taxa_goal=nr_leaves, max_retics=nr_reticulations)
        if nw is not None:
            if nw.nr_of_leaves() == nr_leaves and rets == nr_reticulations:
                # Found a network
                break
    return nw
