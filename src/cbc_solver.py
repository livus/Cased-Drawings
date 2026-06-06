"""
CBC solver implementation for cased drawing optimization using PuLP.
Requires: pulp (free, open-source)
Install: pip install pulp
"""

from itertools import permutations
from typing import Tuple

import gdMetriX as gx
import networkx as nx
from pulp import *

import draw_cd
from draw_cd import EncasedCrossing
from solver_interface import SolverInterface, OptimizationGoal, CasedDrawingModel

logger = logging.getLogger(__name__)


class CBCSolver(SolverInterface):
    """MILP solver using PuLP with CBC backend."""

    def __init__(self):
        """Initialize CBC solver."""
        self.solver = PULP_CBC_CMD(msg=False, timeLimit=1800)

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
        Solve cased drawing problem using CBC solver.

        :param g: Graph component
        :param goal: Optimization goal
        :param model: Cased drawing model
        :param pos: Node positions
        :param time_limit: Time limit in seconds
        :param memory_limit: Memory limit in GB (ignored for CBC)
        :return: Tuple of (encased crossings list, objective value)
        """

        logging.debug("Building CBC model...")

        pos = gx.get_node_positions(g, pos)
        crossings = gx.get_crossings(g, pos)

        if not crossings:
            return [], 0

        big_M = g.number_of_nodes() + g.number_of_edges()

        crossings_per_edge = draw_cd.get_crossings_per_edge_sorted(crossings, pos)
        cr_nr_per_edge = [len(crossings) for crossings in crossings_per_edge.values()]

        # Get index of an edge
        _edge_keys = list(crossings_per_edge.keys())
        edge_index = {_edge_keys[i]: i for i in range(len(_edge_keys))}

        try:
            # Create optimization problem
            prob = LpProblem("CasedDrawing", LpMinimize)

            # Binary variables: c[i,j] = 1 if jth crossing on edge i is a bridge
            c = {}
            for i in range(len(crossings_per_edge)):
                for j in range(cr_nr_per_edge[i]):
                    c[(i, j)] = LpVariable(f"c_{i}_{j}", cat="Binary")

            # Binary variables: s[i,j] = 1 if there's a switch between crossings j and j+1 on edge i
            s = {}
            for i in range(len(crossings_per_edge)):
                for j in range(cr_nr_per_edge[i] - 1):
                    s[(i, j)] = LpVariable(f"s_{i}_{j}", cat="Binary")

            # ===== COMMON CONSTRAINTS =====

            # (1) Exactly one crossing is a bridge at each crossing
            for crossing_index in range(len(crossings)):
                prob += (
                    lpSum(
                        [
                            c[
                                edge_index[edge],
                                crossings_per_edge[edge].index(crossings[crossing_index]),
                            ]
                            for edge in crossings[crossing_index].involved_edges
                        ]
                    )
                    == 1
                ), f"one_bridge_only_{crossing_index}"

            # (2) Linearize switches: s[i,j] = 1 iff c[i,j] != c[i,j+1]
            for i in range(len(cr_nr_per_edge)):
                for j in range(cr_nr_per_edge[i] - 1):
                    prob += c[i, j] - c[i, j + 1] <= s[i, j], f"s_bound_1_{i}_{j}"
                    prob += c[i, j + 1] - c[i, j] <= s[i, j], f"s_bound_2_{i}_{j}"
                    prob += c[i, j] + c[i, j + 1] >= s[i, j], f"s_bound_3_{i}_{j}"
                    prob += (
                        2 - (c[i, j] + c[i, j + 1]) >= s[i, j]
                    ), f"s_bound_4_{i}_{j}"

            # ===== OPTIMIZATION GOALS =====
            if goal == OptimizationGoal.MinTotalSwitches:
                prob += lpSum([s[i, j] for i in range(len(cr_nr_per_edge)) for j in range(cr_nr_per_edge[i] - 1)])

            elif goal == OptimizationGoal.MaxTotalSwitches:
                prob += -lpSum(
                    [s[i, j] for i in range(len(cr_nr_per_edge)) for j in range(cr_nr_per_edge[i] - 1)]
                )

            elif goal == OptimizationGoal.MinMaxSwitches:
                z = LpVariable("z", lowBound=0, cat="Integer")
                for i in range(len(cr_nr_per_edge)):
                    prob += (
                        lpSum([s[i, j] for j in range(cr_nr_per_edge[i] - 1)]) <= z,
                        f"switches_per_edge_{i}",
                    )
                prob += z

            elif goal == OptimizationGoal.MinSwitchEdges:
                b = {}
                for i in range(len(cr_nr_per_edge)):
                    b[i] = LpVariable(f"b_{i}", cat="Binary")
                    prob += (
                        lpSum([s[i, j] for j in range(cr_nr_per_edge[i] - 1)]) <= big_M * b[i],
                        f"edges_with_at_least_one_switch_{i}",
                    )
                prob += lpSum([b[i] for i in range(len(cr_nr_per_edge))])

            # ===== MODEL RESTRICTIONS =====
            if model == CasedDrawingModel.Weaving:
                # No additional constraints
                pass

            elif model == CasedDrawingModel.Stacking:
                # Variables for total order
                o = {}
                for i in range(len(cr_nr_per_edge)):
                    o[i] = LpVariable(f"o_{i}", lowBound=0, upBound=len(cr_nr_per_edge), cat="Integer")

                b2 = {}
                for i in range(len(cr_nr_per_edge)):
                    for j in range(len(cr_nr_per_edge)):
                        if i != j:
                            b2[(i, j)] = LpVariable(f"b2_{i}_{j}", cat="Binary")

                # Total order constraints
                for i in range(len(cr_nr_per_edge)):
                    for j in range(len(cr_nr_per_edge)):
                        if i != j:
                            prob += (
                                o[i] - o[j] - big_M * b2[i, j] <= -1,
                                f"total_order_1_{i}_{j}",
                            )
                            prob += (
                                o[j] - o[i] - big_M * (1 - b2[i, j]) <= -1,
                                f"total_order_2_{i}_{j}",
                            )

                # Crossings must obey total order
                constr_idx = 0
                for i in range(len(cr_nr_per_edge)):
                    for j in range(cr_nr_per_edge[i]):
                        for k in crossings_per_edge[_edge_keys[i]][j].involved_edges:
                            if edge_index[k] != i:
                                prob += (
                                    o[i] - o[edge_index[k]] + big_M * (1 - c[i, j]) >= 0,
                                    f"stacking_model_adherence_{constr_idx}",
                                )
                                constr_idx += 1

            elif model == CasedDrawingModel.Realizable:
                # Height variables for each edge
                d = {}
                l = {}
                for i in range(len(cr_nr_per_edge)):
                    d[i] = LpVariable(f"d_{i}", lowBound=0, upBound=1, cat="Continuous")
                    l[i] = LpVariable(f"l_{i}", lowBound=0, upBound=1, cat="Continuous")

                constr_idx = 0
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

                        prob += (
                            (d[edge_a_idx] + l[edge_a_idx] * edge_a_percentage)
                            - (d[edge_b_idx] + l[edge_b_idx] * edge_b_percentage)
                            + 2 * (1 - c[edge_a_idx, cr_idx])
                            >= 0,
                            f"realizable_model_{constr_idx}",
                        )
                        constr_idx += 1

            # ===== SOLVE =====
            logging.debug("Solving model with CBC...")
            self.solver.timeLimit = time_limit
            prob.solve(self.solver)

            if LpStatus[prob.status] == "Optimal":
                logging.debug(f"Found solution with objective value {int(value(prob.objective))}")

                encased_crossings = []

                for crossing in crossings:
                    encased_crossing = EncasedCrossing.from_crossing(crossing, None)

                    for edge_a_index in range(len(encased_crossing.involved_edges)):
                        edge_a = list(encased_crossing.involved_edges)[edge_a_index]

                        if (
                            value(c[edge_index[edge_a], crossings_per_edge[edge_a].index(crossing)])
                            > 0.5
                        ):
                            encased_crossing.top_edge = edge_a

                    encased_crossings.append(encased_crossing)

                return encased_crossings, int(value(prob.objective))
            else:
                logging.warning(f"No optimal solution found. Status: {LpStatus[prob.status]}")
                return None, sys.maxsize

        except Exception as e:
            logging.error(f"CBC solver error: {e}")
            import traceback
            traceback.print_exc()
            return None, sys.maxsize

