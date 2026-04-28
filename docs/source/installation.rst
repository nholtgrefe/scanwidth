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

Node-scanwidth ILP backend dependencies:

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


Dependencies
============
The following dependencies are required for the base package,
but are automatically installed when installing the package as above.

- `networkx>=3.0.0`
- `numpy>=1.20.0`

for the ILP backend:
- `scipy>=1.9.0`
- `gurobipy>=10.0.0`

for the documentation:
- `sphinx>=7.0.0`
- `pydata-sphinx-theme>=0.15.0`

for the development:
- `pytest>=8.0.0`