"""
scanwidth
=========

`scanwidth` is a Python package for computing the edge- and node-scanwidth of a directed acyclic graph (DAG).

See https://github.com/nholtgrefe/scanwidth for complete documentation.
"""

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth import edge_scanwidth
from scanwidth.extension import Extension
from scanwidth.node_scanwidth import node_scanwidth
from scanwidth.tree_extension import TreeExtension

__all__ = [
    "DAG",
    "Extension",
    "TreeExtension",
    "edge_scanwidth",
    "node_scanwidth",
]

__version__ = "0.2.6"
