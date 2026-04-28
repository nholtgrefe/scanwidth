"""Scanwidth package for computing scanwidth of directed acyclic graphs (DAGs)."""

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

__version__ = "0.2.1"
