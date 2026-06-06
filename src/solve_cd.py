"""
Creates cased drawing under different models and optimization goals.

This module provides solver-agnostic functions for computing cased drawings.
Specific solvers (Gurobi, CBC) can be selected when needed.
"""

from __future__ import annotations

import logging
import sys

logger = logging.getLogger(__name__)

from typing import List, Optional, Dict, Tuple

import gdMetriX
import networkx as nx
import gdMetriX as gx

from draw_cd import EncasedCrossing
from solver_interface import OptimizationGoal, CasedDrawingModel, SolverInterface


def encase_drawing(
    g: nx.Graph,
    goal: OptimizationGoal,
    model: CasedDrawingModel,
    solver: SolverInterface,
    pos: str | Dict | None = None,
    time_limit: int = 1800,
    memory_limit: int = 8,
) -> Tuple[Optional[List[EncasedCrossing]], int]:
    """
    Find a cased drawing for the given embedding.

    :param g: A networkX graph
    :type g: nx.Graph
    :param goal: The optimization goal
    :type goal: OptimizationGoal
    :param model: The model restricting the solution space
    :type model: CasedDrawingModel
    :param solver: MILP solver instance
    :type solver: SolverInterface
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
    total_costs = 0

    components = list(nx.connected_components(planarized_g))

    logger.debug(f"Number of components: {len(components)}")

    for component in components:
        logger.debug(f"Solving component of size {len(component)}")
        component_graph = split_g.subgraph(component)
        casing, cost = _encase_drawing_per_component(
            component_graph,
            goal,
            model,
            solver,
            gx.get_node_positions(component_graph),
            time_limit,
            memory_limit,
        )
        if casing is not None:
            encased_crossings += casing
            total_costs += cost
        else:
            return None, sys.maxsize

    # Replace new node names with original ones again
    for crossing in encased_crossings:
        crossing.involved_edges = [
            (edge[0][0], edge[1][0]) for edge in crossing.involved_edges
        ]
        crossing.top_edge = (crossing.top_edge[0][0], crossing.top_edge[1][0])

    return encased_crossings, total_costs


def _encase_drawing_per_component(
    g: nx.Graph,
    goal: OptimizationGoal,
    model: CasedDrawingModel,
    solver: SolverInterface,
    pos: str | Dict | None,
    time_limit: int,
    memory_limit: int,
) -> Tuple[Optional[List[EncasedCrossing]], int]:
    """
    Solve the cased drawing problem for a single component.

    :param g: Graph component
    :param goal: Optimization goal
    :param model: Cased drawing model
    :param solver: MILP solver instance
    :param pos: Node positions
    :param time_limit: Time limit in seconds
    :param memory_limit: Memory limit in GB
    :return: Tuple of (encased crossings list, objective value)
    """
    return solver.solve(g, goal, model, pos, time_limit, memory_limit)