"""Scanwidth package for computing scanwidth of directed acyclic graphs (DAGs)."""

from scanwidth.dag import DAG
from scanwidth.extension import Extension
from scanwidth.tree_extension import TreeExtension

__all__ = ["DAG", "Extension", "TreeExtension"]
__version__ = "0.1.0"
