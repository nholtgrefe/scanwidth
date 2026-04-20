# scanwidth

Package for computing scanwidth of directed acyclic graphs (DAGs).

## Installation

```bash
pip install scanwidth
```

## Quick Start

```python
from scanwidth import DAG, Extension, TreeExtension

# Load a DAG from a file
dag = DAG("path/to/graph.el")

# Compute scanwidth using various algorithms
sw, extension = dag.optimal_scanwidth()        # Exact algorithm
sw, extension = dag.greedy_heuristic()        # Greedy heuristic
sw, extension = dag.cut_splitting_heuristic()  # Cut-splitting heuristic
sw, extension = dag.simulated_annealing()      # Simulated annealing

# Save the extension
extension.save_file("output.txt")
```

## Algorithms

The package provides several algorithms for computing scanwidth:

- **`optimal_scanwidth()`** - Exact algorithm using dynamic programming
- **`greedy_heuristic()`** - Fast greedy heuristic
- **`cut_splitting_heuristic()`** - Cut-splitting heuristic
- **`simulated_annealing()`** - Simulated annealing metaheuristic

All methods return a tuple `(scanwidth_value, extension_object)`. For specific parameter settings and detailed documentation, refer to the source code.

## Citation

If you use this package in your research, please cite:

**Exact and heuristic computation of the scanwidth of directed acyclic graphs**. *Niels Holtgrefe, Leo van Iersel, and Mark Jones*. Journal of Computer and System Sciences, 160:103802, 2026. doi: [10.1016/j.jcss.2026.103802](https://doi.org/10.1016/j.jcss.2026.103802)

## License

MIT License - see LICENSE file for details.
