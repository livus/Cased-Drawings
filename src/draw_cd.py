from typing import List

import gdMetriX
import gdMetriX as gx
import matplotlib.pyplot as plt
import networkx as nx


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


def draw_edge_casing(
    casing: List[EncasedCrossing], pos=None, tunnel_width=30, tunnel_length=0.5, ax=None
):

    if ax is None:
        ax = plt.gca()

    for crossing in casing:

        if crossing is None or crossing.top_edge is None:
            continue

        edge_vector = gx.Vector.from_point(
            pos[crossing.top_edge[0]]
        ) - gx.Vector.from_point(pos[crossing.top_edge[1]])
        offset = gx.Vector(tunnel_length / 2, 0)
        offset = offset.rotate(edge_vector.rad())

        x = [crossing.pos.x - offset.x, crossing.pos.x + offset.x]
        y = [crossing.pos.y - offset.y, crossing.pos.y + offset.y]

        ax.plot(x, y, linewidth=tunnel_width, color="white", solid_capstyle="butt")

    for crossing in casing:

        if crossing is None or crossing.top_edge is None:
            continue

        edge_vector = gx.Vector.from_point(
            pos[crossing.top_edge[0]]
        ) - gx.Vector.from_point(pos[crossing.top_edge[1]])
        offset = gx.Vector(tunnel_length / 2, 0)
        offset = offset.rotate(edge_vector.rad())

        x_larger = [crossing.pos.x - offset.x*1.1, crossing.pos.x + offset.x*1.1]
        y_larger = [crossing.pos.y - offset.y*1.1, crossing.pos.y + offset.y*1.1]

        ax.plot(x_larger, y_larger, color="black", solid_capstyle="projecting")


def draw_cased_graph(
    g: nx.Graph, casing: List[EncasedCrossing], ax = None, tunnel_width=30, tunnel_length=0.5
) -> None:
    pos = gdMetriX.get_node_positions(g)

    if ax is None:
        ax = plt.gca()

    for edge in g.edges():
        x = [pos[edge[0]][0], pos[edge[1]][0]]
        y = [pos[edge[0]][1], pos[edge[1]][1]]
        ax.plot(x, y, color="black", solid_capstyle="butt")

    draw_edge_casing(casing, pos, tunnel_width, tunnel_length, ax=ax)

    #nx.draw_networkx_nodes(g, pos, ax=ax, node_size=180)
    nx.draw_networkx_nodes(g, pos, ax=ax, node_size=12, node_color="black")
    #nx.draw_networkx_labels(g, pos, ax=ax, font_color="white")
