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

**Paper:** *Exact and Heuristic Computation of the Scanwidth of Directed Acyclic Graphs* by Niels Holtgrefe, Leo van Iersel, and Mark Jones (2024)  
Available at: [arXiv:2403.12734](https://arxiv.org/abs/2403.12734)

## License

MIT License - see LICENSE file for details.
