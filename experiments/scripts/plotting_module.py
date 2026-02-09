import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import time


def plot_graph(graph, figsize):
    """Plots the graph."""
    plt.rcParams['figure.dpi'] = 300

    if figsize != None: plt.figure(1,figsize=figsize, frameon=False) #standard: (6.4, 4.8)
    else: plt.figure(1, frameon=False)
    
    pos = graphviz_layout(graph, prog='dot')
  
    # Draw nodes of tree
    nx.draw_networkx_nodes(graph, pos=pos, node_size=200, node_color='white',
                           edgecolors='black', alpha=1)
    
    # Draw node labels
    nx.draw_networkx_labels(graph, pos=pos, font_size=8)
    
    # Draw graph edges
    nx.draw_networkx_edges(graph, pos=pos, edgelist=list(graph.edges()), 
                           edge_color='black', width=1, arrowsize=11, node_size=200,
                           alpha=0.85, arrowstyle='->')
    plt.show()  


def plot_tree_extension(graph, tree, fancy, offsetting, figsize):
    """Plots the tree-extension with the underlying graph."""
    plt.rcParams['figure.dpi'] = 300

    if figsize != None: plt.figure(1,figsize=figsize, frameon=False) #standard: (6.4, 4.8)
    else: plt.figure(1, frameon=False)
    
    gamma = nx.MultiDiGraph(tree)
    pos = graphviz_layout(gamma, prog='dot')

    # Draw nodes of tree
    nx.draw_networkx_nodes(gamma, pos=pos, node_size=70, node_color='white',edgecolors='black', alpha=1)
    
    # Draw edges of tree
    nx.draw_networkx_edges(gamma, pos=pos, edgelist = list(gamma.edges()),node_size=0,
                           edge_color='lightgray',arrowstyle='-', width=15, alpha=1)
    
    # Draw node labels
    nx.draw_networkx_labels(gamma, pos=pos, font_size=5)
    
    # Draw graph edges
    
    # Parameters
    arrowsize = 9
    width = 1
    edge_color = 'black'
    alpha = 0.7
    offset = offsetting[1]
    init_offset = offsetting[0]
    
    if fancy == False:
        nx.draw_networkx_edges(gamma, pos=pos, edgelist=list(graph.edges()), 
                               edge_color=edge_color, width=width, arrowsize=arrowsize, node_size=60,
                               alpha=alpha, arrowstyle='->')
        plt.show()
        return
    
    all_nodes_names = {v: [v] for v in graph.nodes()}
    all_nodes_pos = {v: [pos[v]] for v in graph.nodes()}

    for u,v in list(graph.edges()):
        route = nx.dijkstra_path(gamma, u, v)
        route_edges = [(route[n],route[n+1]) for n in range(len(route)-1)]
        
        if pos[route[0]] > pos[route[-1]]:
            direction = 'left'
        else:
            direction = 'right'
        
        if len(route_edges) == 1:
            nx.draw_networkx_edges(gamma, pos=pos, edgelist=route_edges, alpha=alpha,
                                    edge_color=edge_color, width=width, arrowsize=arrowsize, 
                                    node_size=60, arrowstyle='->')
            continue
        
        for i, edge in enumerate(route_edges):
            startnode = edge[0]
            endnode = edge[1]
            
            if i == 0:
                new = str(time.time())[-8:-1]
                all_nodes_names[endnode].append(new)
                
                old_pos = all_nodes_pos[endnode][-1]
                
                x_values = [x for (x, y) in all_nodes_pos[endnode]]
                
                if direction == 'left':
                    lefters = [x for x in x_values if x < x_values[0]]
                    extra = init_offset if len(lefters) == 0 else 0
                    new_x = min(x_values) - offset - extra
                elif direction == 'right':
                    righters = [x for x in x_values if x > x_values[0]]
                    extra = init_offset if len(righters) == 0 else 0
                    new_x = max(x_values) + offset + extra
                
                new_pos = (new_x, old_pos[1])
                all_nodes_pos[endnode].append(new_pos)
                pos[new] = new_pos
                
                nx.draw_networkx_nodes(gamma, pos={new: new_pos}, nodelist=[new], node_size=0, alpha=None)
                
                nx.draw_networkx_edges(gamma, pos=pos, edgelist=[(startnode, new)], 
                       edge_color=edge_color, width=width, arrowsize=arrowsize, alpha=alpha,
                       arrowstyle='-', node_size=0)
                
            elif i > 0 and i < len(route_edges)-1:
                new = str(time.time())[-8:-1]
                all_nodes_names[endnode].append(new)
                
                old_pos = all_nodes_pos[endnode][-1]
                
                x_values = [x for (x, y) in all_nodes_pos[endnode]]
                
                if direction == 'left':
                    lefters = [x for x in x_values if x < x_values[0]]
                    extra = init_offset if len(lefters) == 0 else 0
                    new_x = min(x_values) - offset - extra
                elif direction == 'right':
                    righters = [x for x in x_values if x > x_values[0]]
                    extra = init_offset if len(righters) == 0 else 0
                    new_x = max(x_values) + offset + extra
                
                new_pos = (new_x, old_pos[1])
                all_nodes_pos[endnode].append(new_pos)
                pos[new] = new_pos
                
                nx.draw_networkx_nodes(gamma, pos={new: new_pos}, nodelist=[new], node_size=0, alpha=None)
                
                nx.draw_networkx_edges(gamma, pos=pos, edgelist=[(all_nodes_names[startnode][-1], new)], 
                       edge_color=edge_color, width=width, arrowsize=arrowsize, arrowstyle='-',
                       alpha=alpha, node_size=0)
                
            elif i == len(route_edges) - 1:
                starter = all_nodes_names[startnode][-1]

                nx.draw_networkx_edges(gamma, pos=pos, edgelist=[(starter, endnode)], 
                       edge_color=edge_color, width=width, arrowsize=arrowsize, arrowstyle='->', 
                       node_size=[0,60], nodelist=[starter, endnode], alpha=alpha)

    
    plt.show()