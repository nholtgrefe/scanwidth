"""ILP exact solver for node scanwidth with selectable backend."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, Tuple

import networkx as nx
import numpy as np

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth.tree_extension import TreeExtension


class _ScipyBackendSolver:
    """Internal ILP backend using SciPy HiGHS."""

    def solve(
        self,
        dag: DAG,
        time_limit: Optional[float],
        mip_rel_gap: float,
        verbose: bool,
    ) -> SolverResult:
        """Solve node scanwidth with the SciPy backend."""
        try:
            from scipy.optimize import (  # pyright: ignore[reportMissingImports]
                Bounds,
                LinearConstraint,
                milp,
            )
            from scipy.sparse import csr_matrix  # pyright: ignore[reportMissingImports]
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise ImportError(
                "SciPy is required for ILP backend='scipy'. Install with "
                "`pip install scanwidth[scipy]` or `pip install scanwidth[ilp]`.",
            ) from exc

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

        for v in nodes:
            coeffs = {s_idx: 1.0}
            for u in nodes:
                coeffs[z_index[(u, v)]] = coeffs.get(z_index[(u, v)], 0.0) - 1.0
            add_row(coeffs=coeffs, lb=0.0, ub=np.inf)

        coeffs = {x_index[(u, v)]: 1.0 for u in nodes for v in nodes}
        add_row(coeffs=coeffs, lb=float(n - 1), ub=float(n - 1))

        for v in nodes:
            coeffs = {x_index[(u, v)]: 1.0 for u in nodes}
            add_row(coeffs=coeffs, lb=-np.inf, ub=1.0)

        for v in nodes:
            add_row(coeffs={x_index[(v, v)]: 1.0}, lb=0.0, ub=0.0)

        for u, v in edges:
            add_row(coeffs={y_index[(u, v)]: 1.0}, lb=1.0, ub=1.0)

        for u in nodes:
            for v in nodes:
                if u == v:
                    continue
                coeffs = {y_index[(u, v)]: 1.0, y_index[(v, u)]: 1.0}
                add_row(coeffs=coeffs, lb=-np.inf, ub=1.0)

        for v in nodes:
            add_row(coeffs={y_index[(v, v)]: 1.0}, lb=1.0, ub=1.0)

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
            "disp": verbose,
            "mip_rel_gap": mip_rel_gap,
        }
        if time_limit is not None:
            options["time_limit"] = float(time_limit)

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

        extension = TreeExtension(dag=dag, tree=tree).to_extension()
        return SolverResult(
            value=extension.node_scanwidth(),
            extension=extension,
            diagnostics={
                "backend": "scipy",
                "highs_status": int(result.status),
                "highs_message": str(result.message),
                "objective": float(result.fun) if result.fun is not None else None,
            },
        )


class _GurobiBackendSolver:
    """Internal ILP backend using gurobipy/Gurobi."""

    def solve(
        self,
        dag: DAG,
        time_limit: Optional[float],
        mip_rel_gap: float,
        verbose: bool,
    ) -> SolverResult:
        """Solve node scanwidth with the Gurobi backend."""
        try:
            import gurobipy as gp
            from gurobipy import GRB
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise ImportError(
                "gurobipy is required for ILP backend='gurobi'. Install with "
                "`pip install scanwidth[gurobi]` or `pip install scanwidth[ilp]`. "
                "You also need a working Gurobi installation and license.",
            ) from exc

        graph = dag.graph
        nodes: List = list(graph.nodes())
        edges = list(graph.edges())
        n = len(nodes)

        if n == 0:
            empty_extension = TreeExtension(dag=dag, tree=nx.DiGraph()).to_extension()
            return SolverResult(value=0, extension=empty_extension)

        model = gp.Model("node_scanwidth_ilp")
        model.Params.OutputFlag = 1 if verbose else 0
        model.Params.MIPGap = float(mip_rel_gap)
        if time_limit is not None:
            model.Params.TimeLimit = float(time_limit)

        x = model.addVars(nodes, nodes, vtype=GRB.BINARY, name="x")
        y = model.addVars(nodes, nodes, vtype=GRB.BINARY, name="y")
        alpha = model.addVars(
            [(u, v, w) for u in nodes for v in nodes for w in nodes if v != w],
            vtype=GRB.BINARY,
            name="alpha",
        )
        z = model.addVars(nodes, nodes, vtype=GRB.BINARY, name="z")
        s = model.addVar(vtype=GRB.INTEGER, lb=0.0, ub=float(n), name="s")

        model.setObjective(s, GRB.MINIMIZE)

        for v in nodes:
            model.addConstr(s >= gp.quicksum(z[u, v] for u in nodes), name=f"width_{v}")

        model.addConstr(
            gp.quicksum(x[u, v] for u in nodes for v in nodes) == n - 1,
            name="tree_edge_count",
        )

        for v in nodes:
            model.addConstr(
                gp.quicksum(x[u, v] for u in nodes) <= 1,
                name=f"parent_bound_{v}",
            )

        for v in nodes:
            model.addConstr(x[v, v] == 0, name=f"tree_no_self_{v}")

        for u, v in edges:
            model.addConstr(y[u, v] == 1, name=f"reach_edge_{u}_{v}")

        for u in nodes:
            for v in nodes:
                if u == v:
                    continue
                model.addConstr(y[u, v] + y[v, u] <= 1, name=f"antisym_{u}_{v}")

        for v in nodes:
            model.addConstr(y[v, v] == 1, name=f"refl_{v}")

        for u in nodes:
            for v in nodes:
                for w in nodes:
                    if v == w:
                        continue
                    model.addConstr(y[u, w] - alpha[u, v, w] >= 0)
                    model.addConstr(alpha[u, v, w] - y[u, v] <= 0)
                    model.addConstr(alpha[u, v, w] - x[v, w] <= 0)
                    model.addConstr(alpha[u, v, w] - y[u, v] - x[v, w] >= -1)
        for u in nodes:
            for w in nodes:
                if u == w:
                    continue
                model.addConstr(
                    y[u, w] <= gp.quicksum(alpha[u, v, w] for v in nodes if v != w),
                )

        for u, w in edges:
            for v in nodes:
                if u == v:
                    continue
                model.addConstr(z[u, v] - y[u, v] - y[v, w] >= -1)

        model.optimize()
        if model.SolCount == 0:
            raise RuntimeError(
                f"No Gurobi solution found. status={model.Status}, "
                f"runtime={model.Runtime:.3f}s",
            )

        tree = nx.DiGraph()
        tree.add_nodes_from(nodes)
        for u in nodes:
            for v in nodes:
                if u != v and x[u, v].X > 0.5:
                    tree.add_edge(u, v)

        extension = TreeExtension(dag=dag, tree=tree).to_extension()
        return SolverResult(
            value=extension.node_scanwidth(),
            extension=extension,
            diagnostics={
                "backend": "gurobi",
                "gurobi_status": int(model.Status),
                "objective": float(model.ObjVal),
                "runtime_seconds": float(model.Runtime),
            },
        )


@dataclass(frozen=True)
class ILPSolver(Solver):
    """Exact node-scanwidth solver based on an ILP formulation.

    Parameters
    ----------
    backend : Literal["scipy", "gurobi"], optional
        MILP backend used to solve the ILP model.
    time_limit : Optional[float], optional
        Time limit in seconds for the MILP solve.
    mip_rel_gap : float, optional
        Relative MIP gap tolerance for the selected backend.
    verbose : bool, optional
        If True, enable backend output.
    """

    backend: Literal["scipy", "gurobi"] = "scipy"
    time_limit: Optional[float] = None
    mip_rel_gap: float = 0.0
    verbose: bool = False

    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth exactly via ILP."""
        if self.backend == "scipy":
            return _ScipyBackendSolver().solve(
                dag=dag,
                time_limit=self.time_limit,
                mip_rel_gap=self.mip_rel_gap,
                verbose=self.verbose,
            )
        if self.backend == "gurobi":
            return _GurobiBackendSolver().solve(
                dag=dag,
                time_limit=self.time_limit,
                mip_rel_gap=self.mip_rel_gap,
                verbose=self.verbose,
            )
        raise ValueError(f"Unknown ILP backend '{self.backend}'.")


__all__ = ["ILPSolver"]
