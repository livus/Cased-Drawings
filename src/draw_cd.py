from typing import List, Dict

import gdMetriX
import gdMetriX as gx
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from shapely.geometry.linestring import LineString
from shapely.geometry.polygon import Polygon


class EncasedCrossing(gx.crossings.Crossing):
    def __init__(self, pos, involved_edges, top_edge):
        super().__init__(pos, involved_edges)
        self.top_edge = top_edge

    @classmethod
    def from_crossing(cls, crossing: gx.crossings.Crossing, top_edge):
        return cls(crossing.pos, crossing.involved_edges, top_edge)

    def __str__(self):
        return "[{}, edges: {}, top_edge: {}]".format(
            self.pos, sorted(self.involved_edges), self.top_edge
        )


def get_crossings_per_edge_sorted(
    crossings: List[gx.crossings.Crossing], pos
) -> Dict[object, List[gx.crossings.Crossing]]:
    crossings_per_edge = {}
    for crossing in crossings:
        for edge in crossing.involved_edges:
            if edge in crossings_per_edge:
                crossings_per_edge[edge].append(crossing)
            else:
                crossings_per_edge[edge] = [crossing]

    # Order crossings along each edge
    for edge, crossing_edges in crossings_per_edge.items():
        crossings_per_edge[edge] = sorted(
            crossing_edges,
            key=lambda cr: _projection_position(
                pos[edge[0]], pos[edge[1]], (cr.pos.x, cr.pos.y)
            ),
        )

    return crossings_per_edge


def _projection_position(a, b, p):
    a, b, p = np.array(a), np.array(b), np.array(p)
    ab = b - a
    ap = p - a
    return np.dot(ap, ab) / np.dot(ab, ab)


def _draw_polygon(polygon, ax) -> None:
    if polygon.geom_type == "Polygon":
        x, y = polygon.exterior.xy
        ax.fill(x, y, fc="black")
    elif polygon.geom_type == "MultiPolygon":
        for poly in polygon.geoms:
            x, y = poly.exterior.xy
            ax.fill(x, y, fc="black")
    else:
        raise ValueError(f"Invalid geometric shape of type {polygon.geom_type}")


def _convert_edge_to_polygon(edge, pos, edge_width):

    start = pos[edge[0]]
    end = pos[edge[1]]

    line_segment = LineString([start, end])
    return line_segment.buffer(edge_width, cap_style="flat")


def draw_cased_edges(
    g: nx.Graph(),
    casing: List[EncasedCrossing],
    pos=None,
    edge_width: float = 0.005,
    tunnel_width: float = 0.05,
    ax=None,
):

    if ax is None:
        ax = plt.gca()

    crossings_per_edge = get_crossings_per_edge_sorted(casing, pos)

    for edge, crossings in crossings_per_edge.items():
        edge_polygon = _convert_edge_to_polygon(edge, pos, edge_width)

        # Remove tunnel gaps from edge
        for crossing in crossings:
            if crossing.top_edge != edge:
                top_edge = _convert_edge_to_polygon(
                    crossing.top_edge, pos, tunnel_width
                )
                edge_polygon = edge_polygon.difference(top_edge)

        # Draw resulting polygon
        _draw_polygon(edge_polygon, ax)

    # Draw uncrossed edges
    rev = set([tuple(reversed(e)) for e in crossings_per_edge])
    remaining_edges = (set(g.edges()) - set(crossings_per_edge.keys())) - rev
    for edge in remaining_edges:
        edge_polygon = _convert_edge_to_polygon(edge, pos, edge_width)
        _draw_polygon(edge_polygon, ax)


def draw_cased_graph(
    g: nx.Graph,
    casing: List[EncasedCrossing],
    ax=None,
    edge_width: float = 0.005,
    tunnel_width: float = 0.05,
) -> None:
    pos = gdMetriX.normalize_positions(
        g, gdMetriX.get_node_positions(g), preserve_aspect_ratio=False
    )

    if ax is None:
        ax = plt.gca()

    draw_cased_edges(g, casing, pos, edge_width, tunnel_width, ax)
    nx.draw_networkx_nodes(g, pos, ax=ax, node_size=12, node_color="black")
