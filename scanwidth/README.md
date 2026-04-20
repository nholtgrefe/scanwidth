# scanwidth

Package for computing scanwidth of directed acyclic graphs (DAGs).

## Installation

```bash
pip install scanwidth
```

## Quick Start

```python
import networkx as nx
from scanwidth import DAG
from scanwidth.edge_scanwidth import edge_scanwidth

# Load a DAG from a file with NetworkX
graph = nx.read_edgelist(
    "path/to/graph.el",
    create_using=nx.DiGraph,
    nodetype=str,
)
dag = DAG(graph)

# Compute scanwidth using various algorithms
sw, extension = edge_scanwidth(dag, algorithm="xp")                   # Exact XP algorithm
sw, extension = edge_scanwidth(dag, algorithm="greedy")               # Greedy heuristic
sw, extension = edge_scanwidth(dag, algorithm="cut_splitting")        # Cut-splitting heuristic
sw, extension = edge_scanwidth(dag, algorithm="simulated_annealing")  # Simulated annealing
```

## Algorithms

The package provides several algorithms for computing edge scanwidth:

- **`edge_scanwidth(..., algorithm="xp")`** - Exact algorithm using dynamic programming
- **`edge_scanwidth(..., algorithm="greedy")`** - Fast greedy heuristic
- **`edge_scanwidth(..., algorithm="cut_splitting")`** - Cut-splitting heuristic
- **`edge_scanwidth(..., algorithm="simulated_annealing")`** - Simulated annealing metaheuristic
- **`edge_scanwidth(..., algorithm="random")`** - Random extension baseline

The main entry function returns `(scanwidth_value, extension_object)` for all
algorithms. For specific parameter settings and detailed documentation, refer to
the source code.

## Citation

If you use this package in your research, please cite:

**Exact and heuristic computation of the scanwidth of directed acyclic graphs**. *Niels Holtgrefe, Leo van Iersel, and Mark Jones*. Journal of Computer and System Sciences, 160:103802, 2026. doi: [10.1016/j.jcss.2026.103802](https://doi.org/10.1016/j.jcss.2026.103802)

## License

MIT License - see LICENSE file for details.
