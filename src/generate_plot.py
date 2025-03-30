import logging
import matplotlib

logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)-8s %(pathname)s:%(lineno)d - %(asctime)s.%(msecs)03d - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import random

import networkx as nx
from matplotlib import pyplot as plt

import cased_drawings
import draw_cd
from cased_drawings import CasedDrawingModel, OptimizationGoal

random.seed(123)

size = 20

random_graph = nx.Graph()
random_graph.add_node(1, pos=(68.5951, 817.096))
random_graph.add_node(2, pos=(169.928, 817.366))
random_graph.add_node(3, pos=(35.8984, 793.857))
random_graph.add_node(4, pos=(157.489, 787.118))
random_graph.add_node(5, pos=(202.624, 790.074))
random_graph.add_node(6, pos=(43.4646, 764.943))
random_graph.add_node(7, pos=(82.3763, 770.888))
random_graph.add_node(8, pos=(200.192, 778.454))
random_graph.add_node(9, pos=(105.345, 754.134))
random_graph.add_node(10, pos=(182.628, 767.375))
random_graph.add_node(11, pos=(89.402, 736.03))
random_graph.add_node(12, pos=(170.468, 752.243))

random_graph.add_edges_from(
    [
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 7),
        (2, 5),
        (2, 9),
        (2, 12),
        (3, 4),
        (3, 6),
        (4, 5),
        (4, 8),
        (5, 11),
        (5, 8),
        (6, 7),
        (6, 9),
        (6, 11),
        (7, 10),
        (8, 10),
        (9, 11),
        (10, 12),
        (11, 12),
    ]
)


# Enable LaTeX-style fonts
matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "axes.linewidth": 0.25,
    }
)

scale = 2
rows, cols = 3, 4  # Define matrix size
fig, axes = plt.subplots(
    rows, cols, figsize=(cols * scale, rows * scale)
)  # Create subplots


for goal_idx, goal in enumerate([o.value for o in OptimizationGoal]):
    axes[0, goal_idx].set_title(
        f"\\textsc{{{OptimizationGoal(goal).name}}}", fontsize=14
    )

for mod_idx, model in enumerate([m.value for m in CasedDrawingModel]):
    axes[mod_idx, 0].set_ylabel(f"{CasedDrawingModel(model).name}", fontsize=13)

for mod_idx, model in enumerate([m.value for m in CasedDrawingModel]):
    for goal_idx, goal in enumerate([o.value for o in OptimizationGoal]):

        ax = axes[mod_idx, goal_idx]  # Get corresponding axis
        ax.set_xticks([])  # Remove x ticks
        ax.set_yticks([])  # Remove y ticks
        ax.set_aspect("equal")
        print(goal)
        print(model)
        casing = cased_drawings.encase_drawing(
            random_graph, OptimizationGoal(goal), CasedDrawingModel(model)
        )

        print(casing)

        draw_cd.draw_cased_graph(random_graph, casing, ax)


plt.tight_layout()  # Adjust layout

plt.savefig(f"overview.pdf", format="pdf")
plt.show()
