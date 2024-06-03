from pathlib import Path

import streamlit as st
from loguru import logger

from src.backend.simulation import Simulation
from src.backend.simulation_config import SimulationConfig
from src.utils.read_config import read_config

# Here you can find some shitty logic with `st.session_state`
# This was an attempt to keep states of figures while user change fields
# Currently it (and reset button) doesn't work and I haven't figured out why

class Application:
    def __init__(self) -> None:
        logger.info("Creating GUI")
        st.set_page_config(
            page_title="Vesicle Simulation",
            layout="wide",
        )
        self.default_config = read_config(
            Path().resolve() / "config/default_params.yaml"
        )
        self.input_fields = dict()
        self.sim_config = SimulationConfig()
        st.session_state.figs = []
        st.session_state.tabs = []

    def run(self):
        self.col_params, self.col_res = st.columns([2, 3])
        with self.col_params:
            st.header("Simulation Parameters")
            self.read_params()
            self.run_button = st.button("Run simulation", key="update_button")

            """
            self.reset_button = st.button("Reset parameters")
            if self.reset_button:
                self.reset_params()
            """

        with self.col_res:
            st.header("Results")
            if self.run_button:
                self.get_figures()
                self.draw_figures()
            # If we move self.fraw_figures() out of the if statement
            # it should keep figures inact when changing parameters
            # But for now it throws errors

    def get_figures(self):
        #pb = st.progress(0, text="Calculations ...")
        pb = None
        self.sim_config.update(self.input_fields)

        simulation = Simulation(
            sim_config=self.sim_config,
            progress_bar=pb,
        )

        # with st.spinner("Calculation"):
        figures = simulation.run()

        # Downsample data in figures if they are too large
        for fig in figures:
            for trace in fig.data:
                if hasattr(trace, 'x') and len(trace.x) > 1000:  # Example threshold
                    trace.x = trace.x[::10]  # Downsample by taking every 10th point
                if hasattr(trace, 'y') and len(trace.y) > 1000:
                    trace.y = trace.y[::10]

        st.session_state.figs = figures
        st.session_state.tabs = st.tabs(
            [f"fig {i}" for i in range(1, len(st.session_state.figs) + 1)]
        )

    def draw_figures(self):
        for fig, tab in zip(st.session_state.figs, st.session_state.tabs):
            with tab:
                st.plotly_chart(fig)

    def read_params_set(self, params_set: dict) -> dict:
        res_d = {}
        for name, default_value in params_set.items():
            if name not in st.session_state:
                st.session_state.__setitem__(name, default_value)
            kwargs = {"label": name, "value": st.session_state[name]}
            match default_value:
                case bool():  # bool must be the first here since bool value can be casted to number
                    val = st.checkbox(**kwargs)
                case str():
                    val = st.text_input(**kwargs)
                case float():
                    val = st.number_input(format="%0.2e", **kwargs)
                case int():
                    val = st.number_input(**kwargs)
                case _:
                    raise TypeError("Unknown input type")

            res_d[name] = val

        return res_d

    def reset_params(self):
        for params_set in self.default_config.values():
            for name, default_value in params_set.items():
                if name in st.session_state:
                    st.session_state.__setitem__(name, default_value)

    def read_params(self):
        tabs = st.tabs(list(self.default_config.keys()))
        for tab, param_set_name in zip(tabs, self.default_config.keys()):
            with tab:
                self.input_fields[param_set_name] = self.read_params_set(
                    self.default_config[param_set_name]
                )


app = Application()
