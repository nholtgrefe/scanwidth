# scanwidth

Python package for computing edge-scanwidth and node-scanwidth of directed
acyclic graphs (DAGs).

## Installation

Install the base package:

```bash
pip install scanwidth
```

Install optional dependency groups:

```bash
# Development and test dependencies
pip install scanwidth[dev]

# Documentation dependencies
pip install scanwidth[docs]

# SciPy backend for node ILP (algorithm="ilp", backend="scipy")
pip install scanwidth[scipy]

# Gurobi backend for node ILP (algorithm="ilp", backend="gurobi")
pip install scanwidth[gurobi]

# Both ILP backends
pip install scanwidth[ilp]
```

`gurobipy` requires a working Gurobi installation and a valid Gurobi license
(typically commercial, with academic licenses available separately from Gurobi).

Build docs locally:

```bash
sphinx-build -b html docs/source docs/build/html
```

## Quick Start

```python
import networkx as nx
from scanwidth import DAG
from scanwidth.edge_scanwidth import edge_scanwidth
from scanwidth.node_scanwidth import node_scanwidth

graph = nx.read_edgelist(
    "path/to/graph.el",
    create_using=nx.DiGraph,
    nodetype=str,
)
dag = DAG(graph)

esw, ext = edge_scanwidth(dag, algorithm="xp")
nsw, ext = node_scanwidth(dag, algorithm="ilp")
```

## Repository Structure

- `src/scanwidth/`: installable Python package.
- `tests/`: repository test suite (not part of the installed package).
- `experiments/`: experimental materials (not part of the installed package).

For details on experiments, see `experiments/README.md`.

## Citation

If you use this repository in research, please cite:

**Exact and heuristic computation of the scanwidth of directed acyclic
graphs**. *Niels Holtgrefe, Leo van Iersel, and Mark Jones*. Journal of
Computer and System Sciences, 160:103802, 2026. doi:
[10.1016/j.jcss.2026.103802](https://doi.org/10.1016/j.jcss.2026.103802)

## License

MIT License - see `LICENSE`.
