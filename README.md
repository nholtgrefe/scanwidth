# Computing the Scanwidth of DAGs

This repository contains the source code for the Python package `scanwidth` as well as the experimental materials for the accompanying paper:

**Exact and heuristic computation of the scanwidth of directed acyclic graphs**. *Niels Holtgrefe, Leo van Iersel, and Mark Jones*. Journal of Computer and System Sciences, 160:103802, 2026. doi: [10.1016/j.jcss.2026.103802](https://doi.org/10.1016/j.jcss.2026.103802)

If you use the package or any of the other contents in this repository in your research, please cite the paper.

> **Note:** The paper is based on the following MSc Thesis by the first author.
**Computing the Scanwidth of Directed Acyclic Graphs**, *Niels Holtgrefe* (2023)  
Available at: [http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5](http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5)

## Repository Structure

- **`scanwidth/`** - The Python package for computing scanwidth (installable via pip)
  - See `scanwidth/README.md` for package documentation and installation instructions
  
- **`experiments/`** - Experimental materials used for the paper and thesis
  - `scripts/` - Original Python scripts used for experiments
  - `networks/` - Network data files (real and synthetic networks)
  - `results/` - Experimental results
  - See `experiments/README.md` for details

The experiments described in the paper and thesis were performed using the Python scripts in the `experiments/scripts/` folder. Version 0.1.0. of the installable `scanwidth` package is based on these scripts.

---

> The preliminary arXiv version of the paper used a separate repository called `ComputingScanwidth`, which has been deleted. This repository is the most up-to-date version of the code and should be used instead.
