"""
Main script for computing cased drawings with selectable solver.

Usage:
    python cased_drawings.py --solver cbc --goal 1 --model 1
    python cased_drawings.py --solver gurobi --goal 2 --model 3 --file mygraph.json

Available solvers:
    - cbc: Free, open-source CBC solver (via PuLP) [RECOMMENDED] ✓
    - gurobi: Commercial Gurobi solver (requires valid license)
"""

import argparse
import logging
import json
import sys
import matplotlib.pyplot as plt
import networkx as nx

import file_loader
import solve_cd
from solver_interface import OptimizationGoal, CasedDrawingModel
import draw_cd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s"
)
logger = logging.getLogger(__name__)


def get_solver(solver_name: str):
    """Factory function to get solver instance."""
    if solver_name == "gurobi":
        try:
            from gurobi_solver import GurobiSolver
            logger.info("Using Gurobi solver")
            return GurobiSolver()
        except ImportError:
            logger.error("Gurobi solver not available. Install with: pip install gurobipy")
            sys.exit(1)
    elif solver_name == "cbc":
        try:
            from cbc_solver import CBCSolver
            logger.info("Using CBC solver (via PuLP)")
            return CBCSolver()
        except ImportError:
            logger.error("CBC solver not available. Install with: pip install pulp")
            sys.exit(1)
    else:
        logger.error(f"Unknown solver: {solver_name}")
        sys.exit(1)


def load_graph(file_path: str = None) -> nx.Graph:
    """Load graph from JSON file or return demo graph."""
    if file_path:
        try:
            with open(file_path, "r") as f:
                g = json.load(f)
            logger.info(f"Loaded graph from {file_path}")
            return file_loader.load_graph(file_path)
        except FileNotFoundError:
            logger.error(f"File '{file_path}' not found.")
            sys.exit(1)
        except json.JSONDecodeError:
            logger.error(f"Invalid JSON file: {file_path}")
            sys.exit(1)
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
        logger.info("Using demo graph")
        return g


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compute cased drawings with selectable MILP solver",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use CBC solver (free, recommended)
  python cased_drawings.py --solver cbc --goal 1 --model 1
  
  # Use Gurobi solver
  python cased_drawings.py --solver gurobi --goal 2 --model 2
  
  # Load custom graph
  python cased_drawings.py --solver cbc --file mygraph.json
  
  # Save output
  python cased_drawings.py --solver cbc --output result.png
        """
    )

    parser.add_argument(
        "--solver",
        type=str,
        choices=["cbc", "gurobi"],
        default="cbc",
        help="MILP solver to use (default: cbc)"
    )

    parser.add_argument(
        "--file",
        type=str,
        default=None,
        help="Path to JSON file containing graph"
    )

    parser.add_argument(
        "--goal",
        type=int,
        choices=[1, 2, 3, 4],
        default=1,
        help="Optimization goal: 1=MinTotalSwitches (default), 2=MaxTotalSwitches, " +
             "3=MinMaxSwitches, 4=MinSwitchEdges"
    )

    parser.add_argument(
        "--model",
        type=int,
        choices=[1, 2, 3],
        default=1,
        help="Cased drawing model: 1=Weaving (default), 2=Stacking, 3=Realizable"
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save output image"
    )

    parser.add_argument(
        "--time-limit",
        type=int,
        default=1800,
        help="Time limit for solver in seconds (default: 1800)"
    )

    parser.add_argument(
        "--memory-limit",
        type=int,
        default=8,
        help="Memory limit for solver in GB (default: 8, only applicable for Gurobi)"
    )

    args = parser.parse_args()

    # Load solver
    solver = get_solver(args.solver)

    # Load graph
    g = load_graph(args.file)

    # Create goal and model enums
    goal = OptimizationGoal(args.goal)
    model = CasedDrawingModel(args.model)

    logger.info(f"Goal: {goal.name}, Model: {model.name}")
    logger.info("Computing cased drawing...")

    # Solve
    encasing, cost = solve_cd.encase_drawing(
        g,
        goal=goal,
        model=model,
        solver=solver,
        time_limit=args.time_limit,
        memory_limit=args.memory_limit
    )

    if encasing is not None:
        logger.info(f"✓ Solution found with objective value: {cost}")
        logger.info(f"  Number of crossings: {len(encasing)}")

        # Visualize
        draw_cd.draw_cased_graph(g, encasing)

        if args.output:
            plt.savefig(args.output, dpi=150, bbox_inches='tight')
            logger.info(f"Saved output to {args.output}")

        plt.show()
    else:
        logger.error("✗ No solution found within time limit")
        sys.exit(1)


if __name__ == "__main__":
    main()

