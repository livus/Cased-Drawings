import json
import logging

import gdMetriX as gx
import matplotlib.pyplot as plt


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)-8s %(pathname)s:%(lineno)d - %(asctime)s.%(msecs)03d - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import random

import networkx as nx



# File where the execution times are stored
TIME_FILENAME = "execution_times.json"
COST_FILENAME = "solution_size.json"


def load_dic(file_name):
    """Load execution times from file, converting string keys back to tuples."""
    try:
        with open(file_name, "r") as f:
            raw_data = json.load(f)
        return {eval(k): v for k, v in raw_data.items()}  # Convert keys back to tuples
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def save_dic(data, file_name):
    """Save execution times to file, converting tuple keys to strings."""
    with open(file_name, "w") as f:
        json.dump({str(k): v for k, v in data.items()}, f, indent=4)


result_dic = load_dic(TIME_FILENAME)
cost_dic = load_dic(COST_FILENAME)

for size in range(50, 10000, 50):
    print("Size",size)
    for i in range(0, 20):
        random_graph = nx.fast_gnp_random_graph(
            size, random.random() * 0.19 + 0.01, random.randint(1, 10000000)
        )
        spring_layout = nx.spring_layout(random_graph)
        nx.set_node_attributes(random_graph, spring_layout, "pos")

        number_of_crossings = gx.number_of_crossings(random_graph, spring_layout)

        print("Crossings:", number_of_crossings)

        if number_of_crossings <= 150:
            plt.figure(figsize=(10, 8))
            nx.draw(
                random_graph, spring_layout,
                node_size=30,
                with_labels=False,
                edge_color="gray",
                node_color="steelblue",
                alpha=0.7
            )
            plt.axis("equal")
            plt.show()
