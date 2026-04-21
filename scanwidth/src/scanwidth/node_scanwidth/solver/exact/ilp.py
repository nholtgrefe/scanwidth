"""ILP exact solver for node scanwidth using SciPy HiGHS."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth.tree_extension import TreeExtension


@dataclass(frozen=True)
class ILPSolver(Solver):
    """Exact node-scanwidth solver based on an ILP formulation.

    Parameters
    ----------
    time_limit : Optional[float], optional
        Time limit in seconds for the MILP solve.
    mip_rel_gap : float, optional
        Relative MIP gap tolerance for HiGHS.
    verbose : bool, optional
        If True, enable HiGHS output.
    """

    time_limit: Optional[float] = None
    mip_rel_gap: float = 0.0
    verbose: bool = False

    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth exactly via MILP."""
        from scipy.optimize import (  # pyright: ignore[reportMissingImports]
            Bounds,
            LinearConstraint,
            milp,
        )
        from scipy.sparse import csr_matrix  # pyright: ignore[reportMissingImports]

        graph = dag.graph
        nodes: List = list(graph.nodes())
        n = len(nodes)
        edges = list(graph.edges())

        if n == 0:
            empty_extension = TreeExtension(dag=dag, tree=nx.DiGraph()).to_extension()
            return SolverResult(value=0, extension=empty_extension)

        x_index: Dict[Tuple[object, object], int] = {}
        y_index: Dict[Tuple[object, object], int] = {}
        alpha_index: Dict[Tuple[object, object, object], int] = {}
        z_index: Dict[Tuple[object, object], int] = {}

        next_idx = 0
        for u in nodes:
            for v in nodes:
                x_index[(u, v)] = next_idx
                next_idx += 1
        for u in nodes:
            for v in nodes:
                y_index[(u, v)] = next_idx
                next_idx += 1
        for u in nodes:
            for v in nodes:
                for w in nodes:
                    if v != w:
                        alpha_index[(u, v, w)] = next_idx
                        next_idx += 1
        for u in nodes:
            for v in nodes:
                z_index[(u, v)] = next_idx
                next_idx += 1
        s_idx = next_idx
        num_vars = s_idx + 1

        c = np.zeros(num_vars, dtype=float)
        c[s_idx] = 1.0
        integrality = np.ones(num_vars, dtype=int)
        lower = np.zeros(num_vars, dtype=float)
        upper = np.ones(num_vars, dtype=float)
        lower[s_idx] = 0.0
        upper[s_idx] = float(n)
        bounds = Bounds(lb=lower, ub=upper)

        rows: List[int] = []
        cols: List[int] = []
        data: List[float] = []
        lb_rows: List[float] = []
        ub_rows: List[float] = []

        row_idx = 0

        def add_row(coeffs: Dict[int, float], lb: float, ub: float) -> None:
            nonlocal row_idx
            for col_idx, coeff in coeffs.items():
                rows.append(row_idx)
                cols.append(col_idx)
                data.append(coeff)
            lb_rows.append(lb)
            ub_rows.append(ub)
            row_idx += 1

        # 1) s >= sum_u z_uv
        for v in nodes:
            coeffs = {s_idx: 1.0}
            for u in nodes:
                coeffs[z_index[(u, v)]] = coeffs.get(z_index[(u, v)], 0.0) - 1.0
            add_row(coeffs=coeffs, lb=0.0, ub=np.inf)

        # 2) Tree has exactly n-1 edges.
        coeffs = {x_index[(u, v)]: 1.0 for u in nodes for v in nodes}
        add_row(coeffs=coeffs, lb=float(n - 1), ub=float(n - 1))

        # 3) At most one parent per vertex.
        for v in nodes:
            coeffs = {x_index[(u, v)]: 1.0 for u in nodes}
            add_row(coeffs=coeffs, lb=-np.inf, ub=1.0)

        # 4) x[v,v] = 0 to forbid self-loops in the tree.
        for v in nodes:
            add_row(coeffs={x_index[(v, v)]: 1.0}, lb=0.0, ub=0.0)

        # 5) Original edges must be represented as reachable paths: y[u,v] = 1.
        for u, v in edges:
            add_row(coeffs={y_index[(u, v)]: 1.0}, lb=1.0, ub=1.0)

        # 6) Antisymmetry on reachability.
        for u in nodes:
            for v in nodes:
                if u == v:
                    continue
                coeffs = {y_index[(u, v)]: 1.0, y_index[(v, u)]: 1.0}
                add_row(coeffs=coeffs, lb=-np.inf, ub=1.0)

        # 7) Reflexive reachability.
        for v in nodes:
            add_row(coeffs={y_index[(v, v)]: 1.0}, lb=1.0, ub=1.0)

        # 8) Transitivity helper constraints with alpha.
        for u in nodes:
            for v in nodes:
                for w in nodes:
                    if v == w:
                        continue
                    add_row(
                        coeffs={
                            y_index[(u, w)]: 1.0,
                            alpha_index[(u, v, w)]: -1.0,
                        },
                        lb=0.0,
                        ub=np.inf,
                    )
                    add_row(
                        coeffs={
                            alpha_index[(u, v, w)]: 1.0,
                            y_index[(u, v)]: -1.0,
                        },
                        lb=-np.inf,
                        ub=0.0,
                    )
                    add_row(
                        coeffs={
                            alpha_index[(u, v, w)]: 1.0,
                            x_index[(v, w)]: -1.0,
                        },
                        lb=-np.inf,
                        ub=0.0,
                    )
                    add_row(
                        coeffs={
                            alpha_index[(u, v, w)]: 1.0,
                            y_index[(u, v)]: -1.0,
                            x_index[(v, w)]: -1.0,
                        },
                        lb=-1.0,
                        ub=np.inf,
                    )

        for u in nodes:
            for w in nodes:
                if u == w:
                    continue
                coeffs = {y_index[(u, w)]: 1.0}
                for v in nodes:
                    if v == w:
                        continue
                    coeffs[alpha_index[(u, v, w)]] = (
                        coeffs.get(alpha_index[(u, v, w)], 0.0) - 1.0
                    )
                add_row(coeffs=coeffs, lb=-np.inf, ub=0.0)

        # 9) GW membership lower bound for edge boundaries.
        for u, w in edges:
            for v in nodes:
                if u == v:
                    continue
                add_row(
                    coeffs={
                        z_index[(u, v)]: 1.0,
                        y_index[(u, v)]: -1.0,
                        y_index[(v, w)]: -1.0,
                    },
                    lb=-1.0,
                    ub=np.inf,
                )

        constraint_matrix = csr_matrix(
            (np.array(data, dtype=float), (np.array(rows), np.array(cols))),
            shape=(row_idx, num_vars),
        )
        constraints = LinearConstraint(
            A=constraint_matrix,
            lb=np.array(lb_rows, dtype=float),
            ub=np.array(ub_rows, dtype=float),
        )

        options = {
            "presolve": True,
            "disp": self.verbose,
            "mip_rel_gap": self.mip_rel_gap,
        }
        if self.time_limit is not None:
            options["time_limit"] = float(self.time_limit)

        result = milp(
            c=c,
            integrality=integrality,
            bounds=bounds,
            constraints=constraints,
            options=options,
        )

        if result.x is None:
            raise RuntimeError(
                f"HiGHS MILP did not produce a solution. status={result.status}; "
                f"message={result.message}"
            )

        x_sol = result.x
        tree = nx.DiGraph()
        tree.add_nodes_from(nodes)
        for (u, v), idx in x_index.items():
            if u != v and x_sol[idx] > 0.5:
                tree.add_edge(u, v)

        tree_extension = TreeExtension(dag=dag, tree=tree)
        extension = tree_extension.to_extension()
        value = extension.node_scanwidth()
        return SolverResult(
            value=value,
            extension=extension,
            diagnostics={
                "highs_status": int(result.status),
                "highs_message": str(result.message),
                "objective": float(result.fun) if result.fun is not None else None,
            },
        )
