# Computing the Scanwidth of DAGs

This repository contains the source code for the Python package `scanwidth` aimed at computing the scanwidth of directed acyclic graphs (DAGs). The algorithms are to c

**Paper:** *Exact and Heuristic Computation of the Scanwidth of Directed Acyclic Graphs* by Niels Holtgrefe, Leo van Iersel, and Mark Jones (2024)  
Available at: [arXiv:2403.12734](https://arxiv.org/abs/2403.12734)

**Thesis:** *Computing the Scanwidth of Directed Acyclic Graphs* by Niels Holtgrefe (2023)  
Available at: [http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5](http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5)

---

> **Note:** The experiments described in the paper and thesis were performed using the Python scripts in the `experiments/scripts/` folder. These scripts have been refactored and cleaned up to form the basis of the installable `scanwidth` package (version 0.1.0). For details about the experimental materials, see `experiments/README.md`.

---

## Installation

The scanwidth package can be installed via pip:

```bash
pip install scanwidth
```

Or from source:
```bash
pip install -e .
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

> **Note:** The preliminary arXiv version of the paper used a separate repository called `ComputingScanwidth`, which is now deprecated. This repository is the most up-to-date version of the code and should be used instead.
