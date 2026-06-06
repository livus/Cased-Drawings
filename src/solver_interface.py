"""
Abstract interface for MILP solvers used in cased drawing optimization.
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import List, Optional, Tuple, Dict

import networkx as nx

from draw_cd import EncasedCrossing


class OptimizationGoal(Enum):
    """Available optimization goals"""
    MinTotalSwitches = 1
    MaxTotalSwitches = 2
    MinMaxSwitches = 3
    MinSwitchEdges = 4


class CasedDrawingModel(Enum):
    """Available models restricting the solution space"""
    Weaving = 1
    Stacking = 2
    Realizable = 3


class SolverInterface(ABC):
    """Abstract base class for MILP solvers."""

    @abstractmethod
    def solve(
        self,
        g: nx.Graph,
        goal: OptimizationGoal,
        model: CasedDrawingModel,
        pos: Dict,
        time_limit: int = 1800,
        memory_limit: int = 8,
    ) -> Tuple[Optional[List[EncasedCrossing]], int]:
        """
        Solve the cased drawing problem for a single component.

        :param g: Graph component (already split)
        :param goal: Optimization goal
        :param model: Cased drawing model (Weaving, Stacking, Realizable)
        :param pos: Node positions
        :param time_limit: Time limit in seconds
        :param memory_limit: Memory limit in GB
        :return: Tuple of (list of encased crossings, objective value)
        """
        pass

