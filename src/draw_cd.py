from typing import List

import gdMetriX as gx
import matplotlib.pyplot as plt

class EncasedCrossing(gx.crossings.Crossing):
    def __init__(self, pos, involved_edges, top_edge):
        super().__init__(pos, involved_edges)
        self.top_edge = top_edge

    @classmethod
    def from_crossing(cls, crossing: gx.crossings.Crossing, top_edge):
        return cls(crossing.pos, crossing.involved_edges, top_edge)

    def __str__(self):
        return "[{}, edges: {}, top_edge: {}]".format(self.pos, sorted(self.involved_edges), self.top_edge)


def draw_edge_casing(casing: List[EncasedCrossing], pos = None, tunnel_width = 30, tunnel_length = 0.5):

    for crossing in casing:

        if crossing is None or crossing.top_edge is None:
            continue

        edge_vector = gx.Vector.from_point(pos[crossing.top_edge[0]]) - gx.Vector.from_point(pos[crossing.top_edge[1]])
        offset = gx.Vector(tunnel_length/2, 0)
        offset = offset.rotate(edge_vector.rad())

        x = [crossing.pos.x - offset.x, crossing.pos.x + offset.x]
        y = [crossing.pos.y - offset.y, crossing.pos.y + offset.y]

        plt.plot(x, y, linewidth=tunnel_width, color="white", solid_capstyle="butt")
        plt.plot(x,y, color="black", solid_capstyle="projecting")