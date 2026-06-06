import json
import statistics

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline

from solve_cd import OptimizationGoal, CasedDrawingModel

matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "axes.linewidth": 0.25,
    }
)

# File where the execution times are stored
RUNTIME_FILENAME = "execution_times.json"
COST_FILENAME = "solution_size.json"


def load_execution_times(file_name):
    """Load execution times from file, converting string keys back to tuples."""
    try:
        with open(file_name, "r") as f:
            raw_data = json.load(f)
        return {eval(k): v for k, v in raw_data.items()}  # Convert keys back to tuples
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def plot_execution_times(tracker, available_goals, available_models, plot_name):
    num_rows = len(available_models)
    num_cols = len(available_goals)

    fig, axes = plt.subplots(
        num_rows,
        num_cols,
        figsize=(4 * num_cols, 4 * num_rows),
        sharex=True,
        sharey=True,
    )

    if num_rows == 1 and num_cols == 1:
        axes = [[axes]]  # Ensure 2D indexing
    elif num_rows == 1:
        axes = [axes]  # Convert to list of lists
    elif num_cols == 1:
        axes = [[ax] for ax in axes]  # Convert to list of lists

    small_runtime = {}
    large_runtime = {}

    for col_idx, current_goal in enumerate(available_goals):
        for row_idx, current_model in enumerate(available_models):

            ax = axes[row_idx][col_idx]

            # Collect data points
            x_vals = []
            y_vals = []
            for (nm_crossings, goal, model), times in tracker.items():
                if (
                    OptimizationGoal(goal) == current_goal
                    and CasedDrawingModel(model) == current_model
                    and nm_crossings < 150
                    and max(times) < 2000
                ):
                    x_vals.extend(
                        [nm_crossings] * len(times)
                    )  # Repeat integer for each execution time
                    y_vals.extend(times)  # Store execution times

            small_values = [
                y_vals[idx] for idx in range(len(x_vals)) if x_vals[idx] <= 50
            ]
            large_values = [
                y_vals[idx] for idx in range(len(x_vals)) if x_vals[idx] > 50
            ]
            small_runtime[(current_goal, current_model)] = (
                statistics.mean(small_values) if len(small_values) > 0 else 0
            )
            large_runtime[(current_goal, current_model)] = (
                statistics.mean(large_values) if len(large_values) > 0 else 0
            )

            # Scatter plot
            if x_vals:
                ax.set_xlim([0, 150])
                ax.set_ylim([0, 2])
                ax.scatter(x_vals, y_vals, alpha=0.5, color="b", s=0.5)
                # ax.plot(sorted(set(x_vals)), [sum(y for x, y in zip(x_vals, y_vals) if x == xi) / x_vals.count(xi) for xi in sorted(set(x_vals))], "r-", lw=1)  # Mean line

                x_vals_set = sorted(set(x_vals))
                if len(x_vals_set) > 4:
                    y_vals = [
                        statistics.mean(
                            [
                                y_vals[idx]
                                for idx in range(len(x_vals))
                                if x_vals[idx] == x_vals_set[x_val]
                            ]
                        )
                        for x_val in range(len(x_vals_set))
                    ]

                    x_smooth = np.linspace(
                        min(x_vals_set), max(x_vals_set), 300
                    )  # Generate smooth x values
                    spline = make_interp_spline(
                        x_vals_set, y_vals, k=2
                    )  # k=2 for moderate smoothing
                    y_smooth = spline(x_smooth)

                    # Moving average function
                    def moving_average(y, window_size=5):
                        return np.convolve(
                            y, np.ones(window_size) / window_size, mode="valid"
                        )

                    # Apply moving average
                    y_ma = moving_average(y_vals, window_size=8)

                    # ax.plot(x_vals_set, y_vals, "g-", lw=1.5)  # Smoothed average line
                    # ax.plot(x_smooth, y_smooth, "r-", lw=1.5)  # Smoothed average line
                    ax.plot(
                        x_vals_set[: len(y_ma)], y_ma, "r-", lw=1.5
                    )  # Smoothed average line

            ax.set_title(
                f"\\textsc{{{current_goal.name}}}\n{current_model.name}", fontsize=20
            )
            if row_idx == 2:
                ax.set_xlabel("Crossings", fontsize=20)
                ax.tick_params(axis="both", which="major", labelsize=10)
            if col_idx == 0:
                ax.set_ylabel("Execution Time (s)", fontsize=20)
                ax.tick_params(axis="both", which="major", labelsize=10)

    for goal in OptimizationGoal:
        s = f"\\textsc{{{goal.name}}} "
        for model in CasedDrawingModel:
            s += f"& {small_runtime[(goal, model)]:.2f}"
            s += f"& {large_runtime[(goal, model)]:.2f}"
        s += " \\\\"
        print(s)
    print()

    plt.tight_layout()
    plt.savefig(plot_name)


result_dic = load_execution_times(RUNTIME_FILENAME)

total_graphs = 0
for entries in result_dic.items():
    total_graphs += len(entries)
print(f"Total # of graphs: {total_graphs/12}")

print(
    "\\rowcolor{gray!15} \\textbf{Runtime [s]} & Small & Large & Small & Large & Small & Large \\\\"
)
plot_execution_times(
    result_dic, list(OptimizationGoal), list(CasedDrawingModel), "runtime.pdf"
)
print("\\hline")
print("\\hline")
print(
    "\\rowcolor{gray!15} \\textbf{Cost} & Small & Large & Small & Large & Small & Large \\\\"
)
result_dic = load_execution_times(COST_FILENAME)
plot_execution_times(
    result_dic, list(OptimizationGoal), list(CasedDrawingModel), "cost.pdf"
)
