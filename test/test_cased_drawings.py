import random

import pytest
import networkx as nx

import cased_drawings
from cased_drawings import OptimizationGoal, CasedDrawingModel


def test_empty_graph():
    g = nx.Graph()

    crossings, cost = cased_drawings.encase_drawing(
        g, OptimizationGoal.MaxTotalSwitches, CasedDrawingModel.Weaving
    )

    assert len(crossings) == 0
    assert cost == 0


def test_single_crossing():
    edges = [(1, 3), (2, 4)]
    g = nx.Graph()
    g.add_node(1, pos=(0, 0))
    g.add_node(2, pos=(0, 1))
    g.add_node(3, pos=(1, 1))
    g.add_node(4, pos=(1, 0))
    g.add_edges_from(edges)

    crossings, cost = cased_drawings.encase_drawing(
        g, OptimizationGoal.MinTotalSwitches, CasedDrawingModel.Weaving
    )

    assert len(crossings) == 1
    assert crossings[0].top_edge in edges
    assert cost == 0


@pytest.mark.parametrize(
    "goal, model",
    [
        (OptimizationGoal.MinTotalSwitches, CasedDrawingModel.Weaving),
        (OptimizationGoal.MinMaxSwitches, CasedDrawingModel.Weaving),
        (OptimizationGoal.MinSwitchEdges, CasedDrawingModel.Weaving),
        (OptimizationGoal.MinTotalSwitches, CasedDrawingModel.Stacking),
        (OptimizationGoal.MinMaxSwitches, CasedDrawingModel.Stacking),
        (OptimizationGoal.MinSwitchEdges, CasedDrawingModel.Stacking),
        (OptimizationGoal.MinTotalSwitches, CasedDrawingModel.Realizable),
        (OptimizationGoal.MinMaxSwitches, CasedDrawingModel.Realizable),
        (OptimizationGoal.MinSwitchEdges, CasedDrawingModel.Realizable),
    ],
)
def test_three_crossings_on_edge_minimize(goal, model):
    g = nx.Graph()
    g.add_node(1, pos=(0, 0))
    g.add_node(2, pos=(4, 0))
    g.add_node(3, pos=(1, 1))
    g.add_node(4, pos=(1, -1))
    g.add_node(5, pos=(2, 1))
    g.add_node(6, pos=(2, -1))
    g.add_node(7, pos=(3, 1))
    g.add_node(8, pos=(3, -1))
    g.add_edges_from([(1, 2), (3, 4), (5, 6), (7, 8)])

    crossings, cost = cased_drawings.encase_drawing(g, goal, model)

    assert len(crossings) == 3

    if crossings[0].top_edge == (1, 2):
        # Long edge is on top
        assert crossings[0].top_edge == crossings[1].top_edge == crossings[2].top_edge
    else:
        # The three short edges are on top
        assert crossings[1].top_edge != (1, 2)
        assert crossings[2].top_edge != (1, 2)

    assert cost == 0


@pytest.mark.parametrize(
    "model",
    [
        CasedDrawingModel.Weaving,
        CasedDrawingModel.Stacking,
        CasedDrawingModel.Realizable,
    ],
)
def test_three_crossings_on_edge_maximize(model):
    g = nx.Graph()
    g.add_node(1, pos=(0, 0))
    g.add_node(2, pos=(4, 0))
    g.add_node(3, pos=(1, 1))
    g.add_node(4, pos=(1, -1))
    g.add_node(5, pos=(2, 1))
    g.add_node(6, pos=(2, -1))
    g.add_node(7, pos=(3, 1))
    g.add_node(8, pos=(3, -1))
    g.add_edges_from([(1, 2), (3, 4), (5, 6), (7, 8)])

    crossings, cost = cased_drawings.encase_drawing(
        g, OptimizationGoal.MaxTotalSwitches, model
    )

    assert len(crossings) == 3

    assert crossings[0].top_edge != crossings[1].top_edge
    assert crossings[1].top_edge != crossings[2].top_edge
    assert (
        crossings[0].top_edge == (1, 2) and crossings[2].top_edge == (1, 2)
    ) or crossings[1].top_edge == (1, 2)

    assert cost == 2


@pytest.mark.parametrize(
    "model",
    [
        CasedDrawingModel.Weaving,
        CasedDrawingModel.Realizable,
    ],
)
def test_triangle_maximize(model):
    g = nx.Graph()
    g.add_node(1, pos=(0, 0))
    g.add_node(2, pos=(1, 1))
    g.add_node(3, pos=(0, -0.1))
    g.add_node(4, pos=(0.6, 1))
    g.add_node(5, pos=(1, -0.1))
    g.add_node(6, pos=(0.4, 1))
    g.add_edges_from([(1, 2), (3, 4), (5, 6)])

    crossings, cost = cased_drawings.encase_drawing(
        g, OptimizationGoal.MaxTotalSwitches, model
    )

    assert len(crossings) == 3

    assert crossings[0].top_edge != crossings[1].top_edge
    assert crossings[1].top_edge != crossings[2].top_edge
    assert crossings[2].top_edge != crossings[0].top_edge

    assert cost == 3


def test_triangle_maximize_stacking():
    g = nx.Graph()
    g.add_node(1, pos=(0, 0))
    g.add_node(2, pos=(1, 0))
    g.add_node(3, pos=(0, -0.1))
    g.add_node(4, pos=(0.6, 1))
    g.add_node(5, pos=(1, -0.1))
    g.add_node(6, pos=(0.4, 1))
    g.add_edges_from([(1, 2), (3, 4), (5, 6)])

    crossings, cost = cased_drawings.encase_drawing(
        g, OptimizationGoal.MaxTotalSwitches, CasedDrawingModel.Stacking
    )

    assert len(crossings) == 3

    assert (
        (crossings[0].top_edge == crossings[1].top_edge)
        or (crossings[0].top_edge == crossings[2].top_edge)
        or (crossings[1].top_edge == crossings[2].top_edge)
    )

    assert cost == 1


@pytest.mark.parametrize(
    "goal, model, name",
    [
        (
            OptimizationGoal.MinTotalSwitches,
            CasedDrawingModel.Weaving,
            "MinTotalSwitches_Weaving",
        ),
        (
            OptimizationGoal.MinMaxSwitches,
            CasedDrawingModel.Weaving,
            "MinMaxSwitches_Weaving",
        ),
        (
            OptimizationGoal.MinSwitchEdges,
            CasedDrawingModel.Weaving,
            "MinSwitchEdges_Weaving",
        ),
        (
            OptimizationGoal.MaxTotalSwitches,
            CasedDrawingModel.Weaving,
            "MaxTotalSwitches_Weaving",
        ),
        (
            OptimizationGoal.MinTotalSwitches,
            CasedDrawingModel.Stacking,
            "MinTotalSwitches_Stacking",
        ),
        (
            OptimizationGoal.MinMaxSwitches,
            CasedDrawingModel.Stacking,
            "MinMaxSwitches_Stacking",
        ),
        (
            OptimizationGoal.MinSwitchEdges,
            CasedDrawingModel.Stacking,
            "MinSwitchEdges_Stacking",
        ),
        (
            OptimizationGoal.MaxTotalSwitches,
            CasedDrawingModel.Stacking,
            "MaxTotalSwitches_Stacking",
        ),
        (
            OptimizationGoal.MinTotalSwitches,
            CasedDrawingModel.Realizable,
            "MinTotalSwitches_Realizable",
        ),
        (
            OptimizationGoal.MinMaxSwitches,
            CasedDrawingModel.Realizable,
            "MinMaxSwitches_Realizable",
        ),
        (
            OptimizationGoal.MinSwitchEdges,
            CasedDrawingModel.Realizable,
            "MinSwitchEdges_Realizable",
        ),
        (
            OptimizationGoal.MaxTotalSwitches,
            CasedDrawingModel.Realizable,
            "MaxTotalSwitches_Realizable",
        ),
    ],
)
def test_larger_graph(goal, model, name):
    random.seed(3458349)

    size = 20

    random_graph = nx.fast_gnp_random_graph(size, 0.1, random.randint(1, 10000000))
    random_embedding = {
        n: [random.randint(-1000, 1000), random.randint(-1000, 1000)]
        for n in range(0, size + 1)
    }
    nx.set_node_attributes(random_graph, random_embedding, "pos")

    crossings, cost = cased_drawings.encase_drawing(random_graph, goal, model)

    assert crossings is not None
    assert cost >= 0
