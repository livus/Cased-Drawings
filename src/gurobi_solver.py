"""
Gurobi solver implementation for cased drawing optimization.
Requires: gurobipy (commercial license required)
"""

import logging
import sys
from itertools import permutations
from typing import List, Optional, Tuple, Dict

import networkx as nx
import gdMetriX as gx
import gurobipy as gp
from gurobipy import GRB

import draw_cd
from draw_cd import EncasedCrossing
from solver_interface import SolverInterface, OptimizationGoal, CasedDrawingModel

logger = logging.getLogger(__name__)


class GurobiSolver(SolverInterface):
    """MILP solver using Gurobi."""

    def __init__(self):
        """Initialize Gurobi environment."""
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag", 0)
        self.env.start()

    def solve(
        self,
        g: nx.Graph,
        goal: OptimizationGoal,
        model: CasedDrawingModel,
        pos: Dict,
        time_limit: int = 1800,
        memory_limit: int = 8,
    ) -> Tuple[Optional[List[EncasedCrossing]], int]:
        """
        Solve cased drawing problem using Gurobi.

        :param g: Graph component
        :param goal: Optimization goal
        :param model: Cased drawing model
        :param pos: Node positions
        :param time_limit: Time limit in seconds
        :param memory_limit: Memory limit in GB
        :return: Tuple of (encased crossings list, objective value)
        """

        logging.debug("Building Gurobi model...")

        pos = gx.get_node_positions(g, pos)
        crossings = gx.get_crossings(g, pos)

        big_M = g.number_of_nodes() + g.number_of_edges()

        crossings_per_edge = draw_cd.get_crossings_per_edge_sorted(crossings, pos)
        cr_nr_per_edge = [len(crossings) for crossings in crossings_per_edge.values()]

        # Get index of an edge
        _edge_keys = list(crossings_per_edge.keys())
        edge_index = {_edge_keys[i]: i for i in range(len(_edge_keys))}

        with gp.Model("CD", env=self.env) as m:
            try:

                # Binary variables: c[i,j] = 1 if jth crossing on edge i is a bridge
                c = m.addVars(
                    [
                        (i, j)
                        for i in range(len(crossings_per_edge))
                        for j in range(cr_nr_per_edge[i])
                    ],
                    vtype=GRB.BINARY,
                    name="c",
                )

                # Binary variables: s[i,j] = 1 if there's a switch between crossings j and j+1 on edge i
                s = m.addVars(
                    [
                        (i, j)
                        for i in range(len(crossings_per_edge))
                        for j in range(cr_nr_per_edge[i] - 1)
                    ],
                    vtype=GRB.BINARY,
                    name="s",
                )

                # ===== COMMON CONSTRAINTS =====

                # (1) Exactly one crossing is a bridge at each crossing
                m.addConstrs(
                    (
                        (
                            gp.quicksum(
                                c[
                                    edge_index[edge],
                                    crossings_per_edge[edge].index(
                                        crossings[crossing_index]
                                    ),
                                ]
                                for edge in crossings[crossing_index].involved_edges
                            )
                            == 1
                        )
                        for crossing_index in range(len(crossings))
                    ),
                    name="one_bridge_only",
                )

                # (2) Linearize switches: s[i,j] = 1 iff c[i,j] != c[i,j+1]
                m.addConstrs(
                    (
                        (c[i, j] - c[i, j + 1] <= s[i, j])
                        for i in range(len(cr_nr_per_edge))
                        for j in range(cr_nr_per_edge[i] - 1)
                    ),
                    name="s_bound_1",
                )
                m.addConstrs(
                    (
                        (c[i, j + 1] - c[i, j] <= s[i, j])
                        for i in range(len(cr_nr_per_edge))
                        for j in range(cr_nr_per_edge[i] - 1)
                    ),
                    name="s_bound_2",
                )
                m.addConstrs(
                    (
                        (
                            c[i, j] + c[i, j + 1] >= s[i, j]
                            for i in range(len(cr_nr_per_edge))
                            for j in range(cr_nr_per_edge[i] - 1)
                        )
                    )
                )
                m.addConstrs(
                    (
                        (
                            2 - (c[i, j] + c[i, j + 1]) >= s[i, j]
                            for i in range(len(cr_nr_per_edge))
                            for j in range(cr_nr_per_edge[i] - 1)
                        )
                    )
                )

                # ===== OPTIMIZATION GOALS =====
                match goal:
                    case OptimizationGoal.MinTotalSwitches:
                        m.setObjective(
                            gp.quicksum(
                                s[i, j]
                                for i in range(len(cr_nr_per_edge))
                                for j in range(cr_nr_per_edge[i] - 1)
                            ),
                            GRB.MINIMIZE,
                        )

                    case OptimizationGoal.MaxTotalSwitches:
                        m.setObjective(
                            gp.quicksum(
                                s[i, j]
                                for i in range(len(cr_nr_per_edge))
                                for j in range(cr_nr_per_edge[i] - 1)
                            ),
                            GRB.MAXIMIZE,
                        )

                    case OptimizationGoal.MinMaxSwitches:
                        z = m.addVar(vtype=GRB.INTEGER, name="z")
                        m.addConstrs(
                            (
                                gp.quicksum(s[i, j] for j in range(cr_nr_per_edge[i] - 1))
                                <= z
                                for i in range(len(cr_nr_per_edge))
                            ),
                            "switches_per_edge",
                        )
                        m.setObjective(z, GRB.MINIMIZE)

                    case OptimizationGoal.MinSwitchEdges:
                        b = m.addVars(len(cr_nr_per_edge), vtype=GRB.BINARY, name="b")
                        m.addConstrs(
                            (
                                (
                                    gp.quicksum(
                                        s[i, j] for j in range(cr_nr_per_edge[i] - 1)
                                    )
                                    <= big_M * b[i]
                                )
                                for i in range(len(cr_nr_per_edge))
                            ),
                            "edges_with_at_least_one_switch",
                        )
                        m.setObjective(
                            gp.quicksum(b[i] for i in range(len(cr_nr_per_edge))),
                            GRB.MINIMIZE,
                        )

                # ===== MODEL RESTRICTIONS =====
                match model:
                    case CasedDrawingModel.Weaving:
                        # No additional constraints
                        pass

                    case CasedDrawingModel.Stacking:
                        # Variables for total order
                        o = m.addVars(len(cr_nr_per_edge), vtype=GRB.INTEGER, name="o")
                        b2 = m.addVars(
                            [
                                (i, j)
                                for i in range(len(cr_nr_per_edge))
                                for j in range(len(cr_nr_per_edge))
                                if i != j
                            ],
                            vtype=GRB.BINARY,
                            name="b2",
                        )

                        # Total order constraints
                        m.addConstrs(
                            o[i] <= len(cr_nr_per_edge) for i in range(len(cr_nr_per_edge))
                        )

                        m.addConstrs(
                            (
                                o[i] - o[j] - big_M * b2[i, j] <= -1
                                for i in range(len(cr_nr_per_edge))
                                for j in range(len(cr_nr_per_edge))
                                if i != j
                            ),
                            name="total_order_1",
                        )
                        m.addConstrs(
                            (
                                o[j] - o[i] - big_M * (1 - b2[i, j]) <= -1
                                for i in range(len(cr_nr_per_edge))
                                for j in range(len(cr_nr_per_edge))
                                if i != j
                            ),
                            name="total_order_2",
                        )

                        # Crossings must obey total order
                        m.addConstrs(
                            (
                                o[i] - o[edge_index[k]] + big_M * (1 - c[i, j]) >= 0
                                for i in range(len(cr_nr_per_edge))
                                for j in range(cr_nr_per_edge[i])
                                for k in crossings_per_edge[_edge_keys[i]][j].involved_edges
                                if edge_index[k] != i
                            ),
                            name="stacking_model_adherence",
                        )

                    case CasedDrawingModel.Realizable:
                        # Height variables for each edge
                        d = m.addVars(
                            len(cr_nr_per_edge), vtype=GRB.CONTINUOUS, name="d", ub=1
                        )
                        l = m.addVars(
                            len(cr_nr_per_edge), vtype=GRB.CONTINUOUS, name="l", ub=1
                        )

                        for crossing in crossings:
                            for edge_a, edge_b in permutations(crossing.involved_edges, 2):
                                edge_a_idx = edge_index[edge_a]
                                edge_b_idx = edge_index[edge_b]

                                edge_a_percentage = draw_cd.projection_position(
                                    pos[edge_a[0]],
                                    pos[edge_a[1]],
                                    (crossing.pos.x, crossing.pos.y),
                                )
                                edge_b_percentage = draw_cd.projection_position(
                                    pos[edge_b[0]],
                                    pos[edge_b[1]],
                                    (crossing.pos.x, crossing.pos.y),
                                )

                                cr_idx = crossings_per_edge[edge_a].index(crossing)

                                m.addConstr(
                                    (
                                        d[edge_a_idx] + l[edge_a_idx] * edge_a_percentage
                                    )  # Height of edge a at crossing
                                    - (
                                        d[edge_b_idx] + l[edge_b_idx] * edge_b_percentage
                                    )  # Height of edge b at crossing
                                    + 2
                                    * (
                                        1 - c[edge_a_idx, cr_idx]
                                    )  # Only apply constraint if edge_a is a top_edge
                                    >= 0,
                                    name=f"realizable_model_{edge_a_idx}_{edge_b_idx}",
                                )

                # ===== SOLVE =====
                m.write("model.lp")
                m.Params.OutputFlag = 0
                m.Params.TimeLimit = time_limit
                m.Params.SoftMemLimit = memory_limit

                logging.debug("Solving model with Gurobi...")
                m.optimize()

                if m.status == GRB.OPTIMAL:
                    logging.debug(f"Found solution with objective value {int(m.objVal)}")

                    encased_crossings = []

                    for crossing in crossings:
                        encased_crossing = EncasedCrossing.from_crossing(crossing, None)

                        for edge_a_index in range(len(encased_crossing.involved_edges)):
                            edge_a = list(encased_crossing.involved_edges)[edge_a_index]

                            if (
                                c[
                                    edge_index[edge_a],
                                    crossings_per_edge[edge_a].index(crossing),
                                ].X
                                > 0.5
                            ):
                                encased_crossing.top_edge = edge_a

                        encased_crossings.append(encased_crossing)

                    return encased_crossings, int(m.objVal)
                else:
                    logging.warning("No optimal solution found within time limit")
                    return None, sys.maxsize

            except gp.GurobiError as e:
                logging.error(f"Gurobi error ({e.errno}): {e}")
                return None, sys.maxsize

