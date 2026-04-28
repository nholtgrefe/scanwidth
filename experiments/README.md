# Experimental Materials

This directory contains all materials used for the experiments described in the following paper. Please cite the paper if you use these materials.

> Niels Holtgrefe, Leo van Iersel, and Mark Jones. *Exact and heuristic computation of the scanwidth of directed acyclic graphs*. Journal of Computer and System Sciences, 160:103802, 2026. doi: [10.1016/j.jcss.2026.103802](https://doi.org/10.1016/j.jcss.2026.103802)

The paper is based on the following master's thesis:

> Niels Holtgrefe. *Computing the Scanwidth of Directed Acyclic Graphs*.  
> MSc thesis, Delft University of Technology, Delft, The Netherlands, 2023.  
> Available at: [http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5](http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5)

### Notes

1. The scripts in this directory form the basis of `scanwidth` version `0.1.0`. The package has since been improved; for new use, prefer the installable package. These scripts are mainly kept for reproducibility of the original experiments.

2. There was an older repository called `ComputingScanwidth` used for the preliminary arXiv version of the paper; it has been removed, and this repository should be used instead.

## Contents

### **`scripts/`** - Original Python Scripts

This directory contains the Python files used in the project:

- **`compute_scanwidth.py`** - The main file containing all the implemented algorithms. This file contains four classes:
  - `DAG` - Main class for computing scanwidth of directed acyclic graphs
  - `Network` - Specifically for phylogenetic networks (extends DAG functionality)
  - `Extension` - Represents an extension of a graph
  - `TreeExtension` - Represents a tree extension
  
  The class `DAG` (or `Network`) is the important class for the scanwidth-algorithms, and uses `NetworkX` as its graph implementation. The important algorithms described in the paper and thesis are methods of the classes `DAG`/`Network`: `optimal_scanwidth`, `greedy_heuristic`, `cut_splitting_heuristic`, and `simulated_annealing`. All of them return a scanwidth value and an extension object. For the specific parameter-settings of these methods, we refer to the documentation in the file.

- **`ZODS_generator.py`** - Contains the Python code used to generate all synthetic networks. This implements the ZODS-generator according to the birth-hybridization process described in the literature.

- **`plotting_module.py`** - Contains basic plotting-functionality that is imported in the main file. Provides functions for visualizing graphs and tree extensions.

#### How to use the scripts

The main-file `compute_scanwidth.py` is meant to be imported in a python script in the same directory. To illustrate its use, here's an example on how to find the exact value of the scanwidth `sw` of a DAG from a file `in_file`, and how to save the optimal extension in the file `out_file`:

```python
from compute_scanwidth import DAG, Network, TreeExtension, Extension
in_file = r'path\to\dag\in\edge-list\format.txt'

G = DAG(in_file)
sw, extension = G.optimal_scanwidth()

out_file = r'file\path\where\extension\should\be\saved.txt'
extension.save_file(out_file)
```

The different classes also contain some other handy methods, such as the `canonical_tree_extension` method of the `Extension` class, which returns a `TreeExtension` object with the same scanwidth as the extension. For a complete overview of all methods and their uses, we refer to the documentation in the python-file. We note that the methods starting with a `_` are meant for internal use only.

### **`networks/`** - Network Data Files

This directory contains two subfolders with the real and synthetic phylogenetic networks used in the paper and thesis:

- **`real_networks/`** - Real phylogenetic networks. All networks are text-files in the edge-list format (`.el`), where each line contains one arc (e.g. `vertex1 vertex2`). The file names of the real networks are the same as in the original source (see http://phylnet.univ-mlv.fr/).

- **`synthetic_networks/`** - Synthetic networks generated using the ZODS generator. The file names are in the format `a_ZODS_Lxxxx_Ryyyy_zzzzzz.el`, where:
  - `xxxx` denotes the number of leaves
  - `yyyy` denotes the number of reticulations
  - `zzzzzz` denotes the number of the network

All networks are text-files in the edge-list format (`.el`), where each line contains one arc (e.g. `vertex1 vertex2`).

### **`results/`** - Experimental Results

This directory contains Excel files with the complete numerical results of all experiments for each network:

- **`data_real_networks.xlsx`** - Results for real networks
- **`data_synthetic_networks.xlsx`** - Results for synthetic networks

The experiments are explained in Section 6 of the paper / Chapter 6 of the thesis.

