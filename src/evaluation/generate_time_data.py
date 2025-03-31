import json
import logging
import sys
import time

import gdMetriX

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)-8s %(pathname)s:%(lineno)d - %(asctime)s.%(msecs)03d - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import random

import networkx as nx

import cased_drawings
from cased_drawings import CasedDrawingModel, OptimizationGoal


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

for size in range(30, 10000):
    for i in range(0, 20):
        random_graph = nx.fast_gnp_random_graph(
            size, random.random() * 0.19 + 0.01, random.randint(1, 10000000)
        )
        spring_layout = nx.spring_layout(random_graph)
        nx.set_node_attributes(random_graph, spring_layout, "pos")

        nm_crossings = gdMetriX.number_of_crossings(random_graph, spring_layout)
        logging.info(nm_crossings)

        if nm_crossings == 0 or nm_crossings > 150:
            continue

        for mod_idx, model in enumerate([m.value for m in CasedDrawingModel]):
            for goal_idx, goal in enumerate([o.value for o in OptimizationGoal]):

                start = time.perf_counter()
                cost = sys.maxsize
                try:
                    casing, cost = cased_drawings.encase_drawing(
                        random_graph,
                        OptimizationGoal(goal),
                        CasedDrawingModel(model),
                        time_limit=120,
                    )
                except:
                    logging.error("Timeout or invalid model")
                    # continue

                end = time.perf_counter()
                dic_key = (nm_crossings, goal, model)
                if dic_key in result_dic:
                    result_dic[dic_key].append(end - start)
                    cost_dic[dic_key].append(cost)
                else:
                    result_dic[dic_key] = [end - start]
                    cost_dic[dic_key] = [cost]

        save_dic(result_dic, TIME_FILENAME)
        save_dic(cost_dic, COST_FILENAME)
