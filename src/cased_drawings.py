"""
Creates cased drawing under different models and optimization goals using Gurobi.


"""

from __future__ import annotations


import logging

import draw_cd

logger = logging.getLogger(__name__)

import json
from typing import List, Optional, Dict
from enum import Enum
import argparse

import gdMetriX
import gurobipy as gp
import numpy as np
from gurobipy import GRB

import networkx as nx
import gdMetriX as gx
import matplotlib.pyplot as plt

from draw_cd import EncasedCrossing, draw_cased_edges


class OptimizationGoal(Enum):
    """
    Available optimization goals
    """

    MinTotalSwitches = 1
    MaxTotalSwitches = 2
    MinMaxSwitches = 3
    MinSwitchEdges = 4


class CasedDrawingModel(Enum):
    """
    Available models restricting the solution space
    """

    Weaving = 1
    Stacking = 2
    Realizable = 3


def encase_drawing(
    g: nx.Graph,
    goal: OptimizationGoal,
    model: CasedDrawingModel,
    pos: str | Dict | None = None,
    time_limit: int = 1800,
    memory_limit: int = 8,
) -> Optional[List[EncasedCrossing]]:
    """
    Find a cased drawing for the given embedding.

    :param g: A networkX graph
    :type g: nx.Graph
    :param goal: The optimization goal
    :type goal: OptimizationGoal
    :param model: The model restricting the solution space
    :type model: CasedDrawingModel
    :param pos: A list of vertex positions. If pos is None, then the positions are expected as properties in the graph.
    :type pos: object
    :param time_limit: Time limit for the MILP solver in seconds
    :type time_limit: int
    :param memory_limit: Memory limit for the MILP solver in GB
    :type memory_limit: int
    :return: In case the cased drawing could be found in the given time limit, a list of all crossings with the top edge
    at each crossing.
    :rtype: Optional[List[EncasedCrossing]]
    """

    logger.debug("Begin splitting graph into components")

    split_g = nx.Graph()
    pos = gx.get_node_positions(g, pos)

    for edge in g.edges():
        first_node = (edge[0], edge)
        second_node = (edge[1], edge)
        split_g.add_node(first_node, pos=pos[edge[0]])
        split_g.add_node(second_node, pos=pos[edge[1]])
        split_g.add_edge(first_node, second_node)

    planarized_g = split_g.copy()
    gdMetriX.planarize(planarized_g, gx.get_node_positions(planarized_g))

    encased_crossings = []

    components = list(nx.connected_components(planarized_g))

    logger.debug(f"Number of components: {len(components)}")

    for component in components:
        logger.debug(f"Solving component of size {len(component)}")
        component_graph = split_g.subgraph(component)
        encased_crossings += _encase_drawing_per_component(
            component_graph,
            goal,
            model,
            gx.get_node_positions(component_graph),
            time_limit,
            memory_limit,
        )

    # Replace new node names with original ones again
    for crossing in encased_crossings:
        crossing.involved_edges = [
            (edge[0][0], edge[1][0]) for edge in crossing.involved_edges
        ]
        crossing.top_edge = (crossing.top_edge[0][0], crossing.top_edge[1][0])

    return encased_crossings


def _encase_drawing_per_component(
    g: nx.Graph,
    goal: OptimizationGoal,
    model: CasedDrawingModel,
    pos: str | Dict | None,
    time_limit: int,
    memory_limit: int,
) -> Optional[List[EncasedCrossing]]:

    logging.debug("Building model...")

    pos = gx.get_node_positions(g, pos)
    crossings = gx.get_crossings(g, pos)

    big_M = g.number_of_nodes() + g.number_of_edges()

    crossings_per_edge = draw_cd.get_crossings_per_edge_sorted(crossings, pos)
    cr_nr_per_edge = [len(crossings) for crossings in crossings_per_edge.values()]

    # Get index of an edge
    _edge_keys = list(crossings_per_edge.keys())
    edge_index = {_edge_keys[i]: i for i in range(len(_edge_keys))}

    with gp.Env() as env, gp.Model("CD", env=env) as m:

        try:

            # A variable c_ei indicating if the ith crossing on e is a bridge

            c = m.addVars(
                [
                    (i, j)
                    for i in range(len(crossings_per_edge))
                    for j in range(cr_nr_per_edge[i])
                ],
                vtype=GRB.BINARY,
                name="c",
            )

            # A variable s_e indicating the number of switches on an edge e

            s = m.addVars(
                [
                    (i, j)
                    for i in range(len(crossings_per_edge))
                    for j in range(cr_nr_per_edge[i] - 1)
                ],
                vtype=GRB.BINARY,
                name="s",
            )

            # Common constraints

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

            # (2) Sum up crossing configurations to switches

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

            # Optimization goals

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

            # Model restrictions

            match model:
                case CasedDrawingModel.Weaving:
                    # Nothing to do
                    pass

                case CasedDrawingModel.Stacking:

                    # Variables for the total order
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

                    # Constraints for the total order

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

                    # Ensure transitivity (optional, already satisfied by above constraints, might help speed thinks up)
                    # m.addConstrs((b2[i, j] + b2[j, i] == 1 for i in range(len(crossings_per_edge)) for j in
                    #              range(len(crossings_per_edge)) if i != j), name="transitivity_1")
                    # m.addConstrs((b2[i, j] + b2[j, k] - b2[i, k] <= 1 for i in range(len(crossings_per_edge)) for j in
                    #              range(len(crossings_per_edge)) for k in range(len(crossings_per_edge)) if
                    #              i != j and j != k and i != k), name="transitivity_2")

                    # Crossings must oblige to the total order

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

                            edge_a_percentage = draw_cd.projection_position(pos[edge_a[0]], pos[edge_a[1]], (crossing.pos.x, crossing.pos.y))
                            edge_b_percentage = draw_cd.projection_position(pos[edge_b[0]], pos[edge_b[1]], (crossing.pos.x, crossing.pos.y))

                            cr_idx = crossings_per_edge[edge_a].index(crossing)

                            m.addConstr( (d[edge_a_idx] + l[edge_a_idx] * edge_a_percentage)  # Height of edge a at crossing
                                         - (d[edge_b_idx] + l[edge_b_idx] * edge_b_percentage) # Height of edge b at crossing
                                         + 2 * (1 - c[edge_a_idx, cr_idx]) # Only apply constraint if edge_a is a top_edge
                                         >= 0, name=f"realizable_model_{edge_a_idx}_{edge_b_idx}" )

            # Solve

            m.write("model.lp")
            m.Params.OutputFlag = 0
            # m.Params.MIPFocus = 2

            m.Params.TimeLimit = time_limit
            m.Params.SoftMemLimit = memory_limit

            logging.debug("Solving model...")

            m.optimize()

            if m.status == GRB.OPTIMAL:
                logging.debug(f"Found solution with optimal value {int(m.objVal)} ")

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

                return encased_crossings
            else:
                return None

        except gp.GurobiError as e:
            logging.error(f"Error ({e.errno}): {e}")


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser(description="Find cased drawing")
    parser.add_argument("--file", type=str, default="")
    parser.add_argument("--model", type=int, default=1)
    parser.add_argument("--goal", type=int, default=1)
    parser.add_argument("--output", type=str, default="")
    args = parser.parse_args()

    if args.file:
        try:
            with open(args.file, "r") as f:
                g = json.load(f)
        except FileNotFoundError:
            logging.error(f"File '{args.file}' not found.")
        except json.JSONDecodeError:
            logging.error(f"Invalid file structure.")

    else:
        # Demo graph

        g = nx.Graph()

        g.add_node(1, pos=(0, 1))
        g.add_node(2, pos=(1, 2))
        g.add_node(3, pos=(1, 0))
        g.add_node(4, pos=(2, 2))
        g.add_node(5, pos=(2, 0))
        g.add_node(6, pos=(3, 2))
        g.add_node(7, pos=(3, 0))
        g.add_node(8, pos=(4, 1))

        g.add_edges_from([(1, 8), (2, 3), (4, 5), (6, 7), (3, 6)])

    pos = gx.get_node_positions(g)
    logging.debug(args.goal, args.model)
    encasing = encase_drawing(
        g, OptimizationGoal(args.goal), CasedDrawingModel(args.model)
    )
    logging.debug(encasing)

    # Draw cased drawing
    draw_cd.draw_cased_graph(g, encasing)

    if args.output:
        plt.savefig(args.output)

    plt.show()
