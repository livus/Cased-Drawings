"""
File-backed interface: Use any external graph editor to save a file (e.g. GraphML from yEd/Gephi).
The notebook UI reloads the file on demand. Pressing "Solve" will re-load the file from disk,
so you can edit in an external tool and Ctrl+S to update the file; then press Solve.

Supported formats: GraphML (.graphml), GEXF (.gexf), JSON node-link (.json), gpickle (.gpickle), edgelist.

This module provides a compact UI with:
- File path input and a "Load" button
- Simple preview of node/edge counts
- Solver selection (cbc/gurobi)
- "Solve" button which reloads the file and runs the solver

"""

from typing import Optional

import ipywidgets as widgets
import matplotlib.pyplot as plt
import networkx as nx
from IPython.display import clear_output

import solve_cd
import draw_cd
from file_loader import load_graph
from solver_interface import OptimizationGoal, CasedDrawingModel


def get_solver_instance(name: str):
    if name == "gurobi":
        try:
            from gurobi_solver import GurobiSolver
            return GurobiSolver()
        except Exception:
            return None
    else:
        try:
            from cbc_solver import CBCSolver
            return CBCSolver()
        except Exception:
            return None


class FileBackedInterface:
    def __init__(self):
        # Widgets
        self.file_path = widgets.Text(value="", description="File:")
        self.load_btn = widgets.Button(description="Load", button_style="info")
        self.solver_dropdown = widgets.Dropdown(options=[("CBC (free)", "cbc"), ("Gurobi", "gurobi")], value="cbc", description="Solver:")
        self.goal_dropdown = widgets.Dropdown(options=[("Min Total Switches", OptimizationGoal.MinTotalSwitches), ("Max Total Switches", OptimizationGoal.MaxTotalSwitches), ("Min Max Switches", OptimizationGoal.MinMaxSwitches), ("Min Switch Edges", OptimizationGoal.MinSwitchEdges)], value=OptimizationGoal.MinTotalSwitches, description="Goal:")
        self.model_dropdown = widgets.Dropdown(options=[("Weaving", CasedDrawingModel.Weaving), ("Stacking", CasedDrawingModel.Stacking), ("Realizable", CasedDrawingModel.Realizable)], value=CasedDrawingModel.Weaving, description="Model:")
        self.solve_btn = widgets.Button(description="Solve (reload file)", button_style="success")
        self.preview = widgets.HTML()
        self.output = widgets.Output()

        # State
        self.graph: Optional[nx.Graph] = None
        self.solver = None
        self.last_casing = None
        self.last_cost = None
        # drawing params
        self.edge_width_slider = widgets.FloatSlider(value=0.005, min=0.001, max=5, step=0.01, description='Edge width:')
        self.tunnel_width_slider = widgets.FloatSlider(value=0.05, min=0.01, max=50, step=0.1, description='Tunnel width:')
        self.node_size_slider = widgets.FloatSlider(value=12, min=1, max=100, step=0.25, description='Node size:')
        self.show_switches_checkbox = widgets.Checkbox(value=False, description='Show switches')

        # Bind callbacks
        self.load_btn.on_click(self._on_load_clicked)
        self.solve_btn.on_click(self._on_solve_clicked)
        self.solver_dropdown.observe(self._on_solver_changed, names='value')
        # drawing param callbacks (do not trigger solver)
        self.edge_width_slider.observe(self._on_drawing_param_changed, names='value')
        self.tunnel_width_slider.observe(self._on_drawing_param_changed, names='value')
        self.node_size_slider.observe(self._on_drawing_param_changed, names='value')
        self.show_switches_checkbox.observe(self._on_drawing_param_changed, names='value')

        # Initialize solver instance
        self._on_solver_changed({'new': self.solver_dropdown.value})

    def _on_solver_changed(self, change):
        solver_name = change['new'] if isinstance(change, dict) else self.solver_dropdown.value
        self.solver = get_solver_instance(solver_name)
        if self.solver is None:
            self.preview.value = f'<div style="color: red;">Solver "{solver_name}" not available (not installed)</div>'
        else:
            self.preview.value = f'<div style="color: green;">Solver "{solver_name}" ready</div>'

    def _on_load_clicked(self, b):
        path = self.file_path.value.strip()
        with self.output:
            clear_output()
            if not path:
                print("Please enter a file path")
                return
            try:
                g = load_graph(path)
                self.graph = g
                n = g.number_of_nodes()
                m = g.number_of_edges()
                self.preview.value = f"<b>Loaded:</b> {path} <br> Nodes: <b>{n}</b> | Edges: <b>{m}</b>"
                print(f"Loaded graph: nodes={n}, edges={m}")
                # quick draw
                fig, ax = plt.subplots(1,1, figsize=(6,4))
                pos = nx.get_node_attributes(g, 'pos')
                if not pos:
                    pos = None
                nx.draw(g, pos=pos, with_labels=True, ax=ax)
                plt.show()
                # reset last solution when loading a new graph
                self.last_casing = None
                self.last_cost = None
            except Exception as e:
                print(f"Error loading graph: {e}")

    def _on_solve_clicked(self, b):
        path = self.file_path.value.strip()
        with self.output:
            clear_output()
            if not path:
                print("Please enter a file path")
                return
            try:
                g = load_graph(path)
                print(f"Reloaded: {path} (nodes={g.number_of_nodes()}, edges={g.number_of_edges()})")
                # remember loaded graph
                self.graph = g
            except Exception as e:
                print(f"Error loading graph: {e}")
                return

            if self.solver is None:
                print("Solver not initialized or not available")
                return

            # Run solver
            try:
                casing, cost = solve_cd.encase_drawing(g, goal=self.goal_dropdown.value, model=self.model_dropdown.value, solver=self.solver)
                if casing is None:
                    print("No solution found")
                    return
                print(f"Solution found: cost={cost}")
                # store last solution
                self.last_casing = casing
                self.last_cost = cost
                # render using current drawing params
                self._render_solution(g, casing, cost)
            except Exception as e:
                print(f"Error while solving: {e}")

    def create_interface(self):
        # arrange drawing parameter controls
        draw_params = widgets.VBox([self.edge_width_slider, self.tunnel_width_slider, self.node_size_slider, self.show_switches_checkbox])

        left = widgets.VBox([self.file_path, self.load_btn, self.preview])
        middle = widgets.VBox([self.solver_dropdown, self.goal_dropdown, self.model_dropdown, self.solve_btn, widgets.HTML("<hr>"), widgets.HTML("<b>Drawing parameters</b>"), draw_params])
        right = widgets.VBox([self.output])
        return widgets.HBox([left, middle, right])

    def _on_drawing_param_changed(self, change):
        # only redraw the last solution or preview graph without running solver
        if self.last_casing is not None and self.graph is not None:
            # redraw solution
            self._render_solution(self.graph, self.last_casing, self.last_cost)
        elif self.graph is not None:
            # redraw preview with new edge width (fallback to simple draw)
            with self.output:
                clear_output()
                fig, ax = plt.subplots(1,1, figsize=(6,4))
                pos = nx.get_node_attributes(self.graph, 'pos')
                if not pos:
                    pos = None
                nx.draw(self.graph, pos=pos, with_labels=True, ax=ax)
                plt.show()

    def _render_solution(self, g, casing, cost):
        # Draw the cased drawing with current drawing parameters and optionally overlay switches
        with self.output:
            clear_output()
            fig, ax = plt.subplots(1,1, figsize=(8,6))
            edge_w = float(self.edge_width_slider.value)
            tunnel_w = float(self.tunnel_width_slider.value)
            node_w = float(self.node_size_slider.value)
            include_switches = self.show_switches_checkbox.value
            draw_cd.draw_cased_graph(g, casing, ax=ax, edge_width=edge_w, tunnel_width=tunnel_w, node_size=node_w, draw_switches=include_switches)

            ax.set_title(f"Cased Drawing - Cost: {cost}")
            plt.tight_layout()
            plt.show()


def create_file_backed_interface():
    ui = FileBackedInterface()
    return ui.create_interface()

