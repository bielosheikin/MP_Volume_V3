from typing import Optional

from streamlit.delta_generator import DeltaGenerator

from src.backend.simulation_config import SimulationConfig
from src.trash import (
    simulation_tools,
)
from src.backend.build_graphs import figure_plotting, plot_dependency


class SimulationResult:
    def __init__(self) -> None:
        pass


class Simulation:
    def __init__(
        self,
        *,
        sim_config: SimulationConfig,
        progress_bar: Optional[DeltaGenerator] = None,
    ):
        self.config = sim_config
        self.pb = progress_bar

    def update_config(self, config: SimulationConfig):
        self.config = config

    def run(self):
        # Do something
        if not self.config:
            raise ValueError("No config provided")
        res_calculate: SimulationResult = self.calculate(self.config)
        result = self.illustrate(res_calculate)

        return result

    def calculate(self, config: SimulationConfig) -> SimulationResult:
        result = simulation_tools.run_simulation(
            config.additional_parameters["internal_ions_amounts"],
            config.parameters,
            progress_bar=self.pb,
        )
        return result

    def illustrate(self, result: SimulationResult):
        # figure1 = display_simulation_results.figure_plottting(
        #     result, self.config.parameters["T"], self.config.parameters["dt"]
        # )
        figure1 = figure_plotting(
            result, self.config.parameters["T"], self.config.parameters["dt"]
        )
        figure2 = plot_dependency()
        return figure1, figure2
