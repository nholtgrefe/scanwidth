Quickstart
==========

`scanwidth` revolves around three core classes:

- :class:`scanwidth.DAG`: validated DAG wrapper around a `networkx.DiGraph`.
- :class:`scanwidth.Extension`: a topological ordering representation used by
  scanwidth computations.
- :class:`scanwidth.TreeExtension`: a rooted tree-extension representation that
  can be converted to an :class:`scanwidth.Extension`.

Most users start from two public API functions:

- :func:`scanwidth.edge_scanwidth` for edge-scanwidth.
- :func:`scanwidth.node_scanwidth` for node-scanwidth.

Below is a minimal end-to-end example that runs both functions on a small DAG
with two sources that merge into one sink.

.. code-block:: python

   import networkx as nx
   from scanwidth import DAG, edge_scanwidth, node_scanwidth

   graph = nx.DiGraph([(1, 3), (2, 3)])
   dag = DAG(graph)

   esw, edge_extension = edge_scanwidth(dag, algorithm="xp")
   nsw, node_extension = node_scanwidth(dag, algorithm="ilp", backend="scipy")

   print(esw, edge_extension)
   print(nsw, node_extension)

What happens here:

1. A `networkx.DiGraph` is wrapped as a :class:`scanwidth.DAG`.
2. :func:`scanwidth.edge_scanwidth` computes edge-scanwidth with the exact XP
   algorithm.
3. :func:`scanwidth.node_scanwidth` computes node-scanwidth via ILP using the
   SciPy backend.
4. Both functions return `(value, extension)`.
5. The returned extension object can also be inspected for its ordering and
   converted to a tree extension when needed.

Reductions
----------

Both :func:`scanwidth.edge_scanwidth` and :func:`scanwidth.node_scanwidth`
support reduction rules by default (``reduce=True``). Reduction often speeds up exact
solvers and usually helps heuristics as well. Parallelization is also supported
by passing `reducer_config=ReducerConfig(parallel_sblocks=True)`.

For complete signatures, keyword arguments, solver classes, reducer
configurations, and class methods, see :doc:`api/index`.

For backend-specific installation details, see :doc:`installation`.

