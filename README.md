[![PyPI](https://img.shields.io/pypi/v/scanwidth)](https://pypi.org/project/scanwidth/)
[![License](https://img.shields.io/github/license/nholtgrefe/scanwidth)](https://github.com/nholtgrefe/scanwidth/blob/main/LICENSE)
[![Docs](https://img.shields.io/badge/docs-stable-blue)](https://nholtgrefe.github.io/scanwidth/)
[![JCSS DOI](https://img.shields.io/badge/JCSS-10.1016%2Fj.jcss.2026.103802-blue)](https://doi.org/10.1016/j.jcss.2026.103802)

# scanwidth

`scanwidth` is a Python package for computing edge-scanwidth and node-scanwidth
of directed acyclic graphs (DAGs). It provides exact and heuristic algorithms,
plus reduction pipelines that make practical computation on larger instances
more tractable.

## Key Features

- **Tree-extension and Extension classess**: classes that support
(tree-)extensions of DAGs, provide computations for scanwidth-bags and 
converting to canonical tree-extenions
- **Exact and heuristic solvers**: XP, brute-force, partition-based exact
  methods (edge), ILP backend selection (node), and multiple heuristics for edge- and node-scanwidth.
- **Reduction framework**: configurable `ReducerConfig`/`Reducer` pipelines for
  edge and node scanwidth, including optional parallel s-block solving.

## Installation

Install the base package:

```bash
pip install scanwidth
```

Install optional ILP dependencies:

```bash
# SciPy backend for node_scanwidth(..., algorithm="ilp", backend="scipy")
pip install scanwidth[scipy]

# Gurobi backend for node_scanwidth(..., algorithm="ilp", backend="gurobi")
pip install scanwidth[gurobi]

# Both ILP backends
pip install scanwidth[ilp]
```

`gurobipy` requires a working Gurobi installation and a valid Gurobi license
(commercial or academic, depending on your setup).

Version requirements for optional ILP backends:

- `scipy>=1.9.0`
- `gurobipy>=10.0.0`

## Documentation

For installation instructions, quickstart examples, and full API reference, see
the **[scanwidth docs](https://nholtgrefe.github.io/scanwidth/)**.

## Citation

If you use `scanwidth` in your research, please cite the corresponding paper:

> Niels Holtgrefe, Leo van Iersel, and Mark Jones. *Exact and heuristic computation of the scanwidth of directed acyclic graphs*. Journal of Computer and System Sciences, 160:103802, 2026. doi: [10.1016/j.jcss.2026.103802](https://doi.org/10.1016/j.jcss.2026.103802)

To view the experimental materials of the paper, go to the folder `experiments`.
