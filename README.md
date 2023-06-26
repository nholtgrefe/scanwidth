# Computing the Scanwidth of DAGs

Repository corresponding to the MSc Thesis "Computing the Scanwidth of Directed Acyclic Graphs" by Niels Holtgrefe, available at https://repository.tudelft.nl/.

Citation in bibtex format:

@masterthesis{holtgrefe2023scanwidth, 
title={Computing the Scanwidth of Directed Acyclic Graphs}, author={Niels Holtgrefe}, year={2023}, month={7}, address={Delft, The Netherlands}, note={Available at \url{https://repository.tudelft.nl/}}, school={Delft University of Technology}, type= {Master's thesis}
}



## Structure
The repository contains three folders:
* 'code', which contains the Python files used in the project. The file 'ZODS_generator.py' contains the python-code used to generate all synthetic networks, while the file 'plotting_module.py' contains some basic plotting-functionality that is imported in the main file. The main file is 'scanwidth.py' which contains all the implemented algorithms. It contains a class Network (and the more general clas DAG) with various methods to call the algorithms. We refer to the documentation in the file on how to run it.
* 'networks', containing two subfolders with the real and synthetic phylogenetic networks used in the thesis. All networks are in edge-list format, where each line of the text file contains one arc (e.g. 'vertex1 vertex2'). The file names of the real networks are the same as in the original source (http://phylnet.univ-mlv.fr/). The file names of the synthetic networks are in the format 'a_ZODS_Lxxxx_Ryyyy_zzzzzz.el', where 'x' denotes the number of leaves, 'y' the number of reticulations, and 'z' the number of the network.
* 'experimental_results', which contains an xlsx-file for the real networks, and one for the synthetic networks, containing the complete numerical results of all experiments for each network.
