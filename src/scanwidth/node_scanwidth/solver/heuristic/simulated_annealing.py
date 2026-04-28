"""Simulated annealing solver for node scanwidth."""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Union

import networkx as nx
import numpy as np

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth.node_scanwidth.solver.heuristic.greedy import GreedySolver
from scanwidth.node_scanwidth.solver.heuristic.random import RandomSolver
from scanwidth.node_scanwidth.solver._utils import node_bag_size
from scanwidth.tree_extension import TreeExtension


@dataclass(frozen=True)
class SimulatedAnnealingSolver(Solver):
    """Simulated-annealing node-scanwidth solver."""

    max_iter: int = 100
    p_in: float = 0.9
    p_stop: float = 0.01
    init_ext: Union[str, Extension] = "greedy"
    verbose: bool = False
    seed: int = 42

    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth using simulated annealing."""
        graph = dag.graph
        rng = np.random.RandomState(self.seed)

        extension = self._initial_extension(dag)
        tree_extension = extension.canonical_tree_extension()
        tree = tree_extension.tree

        sw_values: Dict = {}
        for v in graph.nodes:
            sw_values[v] = node_bag_size(graph, nx.descendants(tree, v) | {v})

        iteration_counter = 0
        best_tree = tree.copy()
        best_sw_values = sw_values.copy()
        best_sw = max(best_sw_values.values())
        epoch_length = len(graph.nodes())

        random_scanwidths: List[int] = []
        for _ in range(100):
            sample = _random_neighbour_scanwidth(
                graph=graph, tree=tree, sw_values=sw_values, rng=rng
            )
            if sample is not None:
                random_scanwidths.append(max(sample[0].values()))

        deltas = [val - best_sw for val in random_scanwidths if val - best_sw > 0]
        delta_f = 1 if not deltas else sum(deltas) / len(deltas)

        init_temp = -delta_f / np.log(self.p_in)
        stop_temp = -delta_f / np.log(self.p_stop)
        temp = init_temp
        vals: List[int] = [max(sw_values.values())]

        while temp >= stop_temp:
            for _ in range(epoch_length):
                res = _random_neighbour_scanwidth(
                    graph=graph, tree=tree, sw_values=sw_values, rng=rng
                )
                if res is None:
                    if self.verbose:
                        print("Only one possible tree extension")
                    return SolverResult(
                        value=max(sw_values.values()),
                        extension=extension,
                        history=[max(sw_values.values()) for _ in range(self.max_iter + 1)],
                    )

                new_sw_values, vertex, parent, connected_vertices = res
                diff = max(new_sw_values.values()) - max(sw_values.values())
                metropolis = np.exp(-diff / temp)

                if diff < 0 or rng.random() < metropolis:
                    sw_values = new_sw_values

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

                    if max(sw_values.values()) <= best_sw:
                        best_tree = tree.copy()
                        best_sw_values = sw_values.copy()
                        best_sw = max(best_sw_values.values())

                iteration_counter += 1

                if self.verbose and iteration_counter % 100 == 0:
                    sys.stdout.write(
                        "\r"
                        + f"current scanwidth: {max(sw_values.values())}; "
                        + f"iteration {iteration_counter}; "
                        + f"best scanwidth: {best_sw}"
                    )
                    time.sleep(1e-40)

            vals.append(max(sw_values.values()))
            alpha = (stop_temp / init_temp) ** (1 / (self.max_iter - 1))
            temp = alpha * temp

        best_tree_extension = TreeExtension(dag, best_tree)
        best_extension = best_tree_extension.to_extension()
        return SolverResult(value=best_sw, extension=best_extension, history=vals)

    def _initial_extension(self, dag: DAG) -> Extension:
        """Return the initial extension for annealing."""
        if isinstance(self.init_ext, Extension):
            return self.init_ext
        if self.init_ext == "greedy":
            return GreedySolver().solve(dag).extension
        if self.init_ext == "random":
            return RandomSolver(seed=self.seed).solve(dag).extension
        raise ValueError(
            f"init_ext must be an Extension object or one of "
            f"{{'greedy', 'random'}}, got {self.init_ext}"
        )


def _random_neighbour_scanwidth(
    graph: nx.DiGraph,
    tree: nx.DiGraph,
    sw_values: Dict,
    rng: np.random.RandomState,
) -> Optional[Tuple[Dict, object, object, Set]]:
    """Return node-scanwidth data for a random neighboring tree extension."""
    possible_choices = []
    for vertex in graph.nodes:
        pred = list(tree.predecessors(vertex))
        if len(pred) == 0:
            continue
        if (pred[0], vertex) in graph.edges:
            continue
        possible_choices.append(vertex)

    if len(possible_choices) == 0:
        return None

    vertex = rng.choice(possible_choices)
    parent = list(tree.predecessors(vertex))[0]

    new_sw_values = sw_values.copy()
    new_sw_values[vertex] = sw_values[parent]

    connected_vertices = {parent}
    for child in tree.successors(parent):
        if child != vertex:
            connected_vertices = connected_vertices | nx.descendants(tree, child) | {child}

    for child in tree.successors(vertex):
        sinkset = nx.descendants(tree, child) | {child}
        connected = False
        for node in sinkset:
            if (parent, node) in graph.edges():
                connected = True
                break
        if connected:
            connected_vertices = connected_vertices | nx.descendants(tree, child) | {child}

    new_sw_values[parent] = node_bag_size(graph, connected_vertices)
    return new_sw_values, vertex, parent, connected_vertices
