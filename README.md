# Computing the Scanwidth of DAGs

This repository contains the source code for the Python package `scanwidth` as well as the experimental materials for the accompanying paper:

**Paper:** *Exact and Heuristic Computation of the Scanwidth of Directed Acyclic Graphs* by Niels Holtgrefe, Leo van Iersel, and Mark Jones (2024)  
Available at: [arXiv:2403.12734](https://arxiv.org/abs/2403.12734)

The paper is based on the following MSc Thesis by the first author.

**Thesis:** *Computing the Scanwidth of Directed Acyclic Graphs* by Niels Holtgrefe (2023)  
Available at: [http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5](http://resolver.tudelft.nl/uuid:9c82fd2a-5841-4aac-8e40-d4d22542cdf5)

## Repository Structure

- **`scanwidth/`** - The Python package for computing scanwidth (installable via pip)
  - See `scanwidth/README.md` for package documentation and installation instructions
  
- **`experiments/`** - Experimental materials used for the paper and thesis
  - `scripts/` - Original Python scripts used for experiments
  - `networks/` - Network data files (real and synthetic networks)
  - `results/` - Experimental results
  - See `experiments/README.md` for details

---

> **Note:** The experiments described in the paper and thesis were performed using the Python scripts in the `experiments/scripts/` folder. Version 0.1.0. of the installable `scanwidth` package is based on these scripts.

---

## Citation

If you use the package or any of the other contents in this repository in your research, please cite:

*Exact and Heuristic Computation of the Scanwidth of Directed Acyclic Graphs* by Niels Holtgrefe, Leo van Iersel, and Mark Jones (2024)  
Available at: [arXiv:2403.12734](https://arxiv.org/abs/2403.12734)

---

> **Note:** The preliminary arXiv version of the paper used a separate repository called `ComputingScanwidth`, which is now deprecated. This repository is the most up-to-date version of the code and should be used instead.
