Installation
============

Base package
------------

Install the core package:

.. code-block:: bash

   pip install scanwidth

Optional dependency groups
--------------------------

Development and testing tools:

.. code-block:: bash

   pip install scanwidth[dev]

Documentation dependencies:

.. code-block:: bash

   pip install scanwidth[docs]

Node-scanwidthILP backend dependencies:

.. code-block:: bash

   # SciPy backend for node_scanwidth(..., algorithm="ilp", backend="scipy")
   pip install scanwidth[scipy]

   # Gurobi backend for node_scanwidth(..., algorithm="ilp", backend="gurobi")
   pip install scanwidth[gurobi]

   # Both ILP backends
   pip install scanwidth[ilp]

- `gurobipy` requires a working Gurobi installation.
- A valid Gurobi license is required (commercial or academic, depending on use).


Development dependencies:
-------------------------

.. code-block:: bash

   pip install scanwidth[dev]
