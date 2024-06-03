"""
Nikita: I only rewrote matplotlib graph building into plotly.
Didn't touch any data composing logic
"""

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
from src.trash import dep_functions as dep_funct

FIG_HEIGHT = 1200
FIG_WIDTH = 1200


def figure_plotting(result, T, dt) -> go.Figure:
    t_axis2 = np.arange(0, T - dt, dt)
    t_axis = np.arange(0, T, dt)
    titles = [
        "pH",
        "Volume",
        "Membrane potential",
        "CL<sup>--</sup> flux through ASOR",
        "Cl<sup>--</sup> flux through CLC$",
        "H<sup>+</sup> flux through CLC",
        "Na<sup>+</sup> flux through TPC",
        "H<sup>+</sup> flux through V-ATPase",
        "H<sup>+</sup> flux through H-leak",
    ]
    fig = make_subplots(
        rows=3,
        cols=3,
        subplot_titles=titles,
        # shared_xaxes=True,
    )

    # pH
    position = {"row": 1, "col": 1}
    fig.add_trace(
        go.Scatter(
            x=t_axis,
            y=result["other_variables"]["vesicle_parameters"]["pH"],
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="pH", range=[4.5, 8.5], **position)

    # Volume
    position = {"row": 1, "col": 2}
    fig.add_trace(
        go.Scatter(
            x=t_axis,
            y=result["other_variables"]["vesicle_parameters"]["V"] * 1e18,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="µm<sup>3</sup>", range=[1, 10], **position)

    # Membrane potential
    position = {"row": 1, "col": 3}
    fig.add_trace(
        go.Scatter(
            x=t_axis,
            y=result["other_variables"]["vesicle_parameters"]["U"] * 1000,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="mV", range=[-90, 60], **position)

    # Cl$^{-}$ flux through ASOR
    position = {"row": 2, "col": 1}
    fig.add_trace(
        go.Scatter(
            x=t_axis2,
            y=result["fluxes"]["Cl"]["ASOR"] * 1e18,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="mol*s<sup>-1</sup>, 10<sup>-18</sup>", **position)

    # Cl$^{-}$ flux throug CLC
    position = {"row": 2, "col": 2}
    fig.add_trace(
        go.Scatter(
            x=t_axis2,
            y=result["fluxes"]["Cl"]["CLC"] * 1e19,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="mol*s<sup>-1</sup>, 10<sup>-19</sup>", **position)

    # H$^{+}$ flux through CLC
    position = {"row": 2, "col": 3}
    fig.add_trace(
        go.Scatter(
            x=t_axis2,
            y=result["fluxes"]["H"]["CLC"] * 1e19,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="mol*s<sup>-1</sup>, 10<sup>-19</sup>", **position)

    # Na$^{+}$ flux through TPC
    position = {"row": 3, "col": 1}
    fig.add_trace(
        go.Scatter(
            x=t_axis2,
            y=result["fluxes"]["Na"]["TPC"] * 1e18,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="mol*s<sup>-1</sup>, 10<sup>-18</sup>", **position)

    # H$^{+}$ flux through V-ATPase
    position = {"row": 3, "col": 2}
    fig.add_trace(
        go.Scatter(
            x=t_axis2,
            y=result["fluxes"]["H"]["vATPase"] * 1e20,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="mol*s<sup>-1</sup>, 10<sup>-20</sup>", **position)

    # H$^{+}$ flux through H-leak
    position = {"row": 3, "col": 3}
    fig.add_trace(
        go.Scatter(
            x=t_axis2,
            y=result["fluxes"]["H"]["H_leak"] * 1e20,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="Time, s", range=[0, T], **position)
    fig.update_yaxes(title_text="mol*s<sup>-1</sup>, 10<sup>-20</sup>", **position)

    fig.update_layout(height=FIG_HEIGHT, width=FIG_WIDTH)

    return fig


def plot_dependency() -> go.Figure:
    pH_step = 0.1
    pH_start = 1
    pH_end = 12
    voltage_step = 5 * 1e-3
    voltage_start = -100 * 1e-3
    voltage_end = 100 * 1e-3

    pH_axis = np.arange(pH_start, pH_end, pH_step)
    voltage_axis = np.arange(voltage_start, voltage_end, voltage_step)

    pH_values_ASOR_wt = np.zeros(len(pH_axis))
    pH_values_ASOR_mutant = np.zeros(len(pH_axis))
    U_values_ASOR = np.zeros(len(voltage_axis))
    pH_values_ClC = np.zeros(len(pH_axis))
    U_values_ClC = np.zeros(len(voltage_axis))

    pH = pH_start

    for i in range(len(pH_axis)):
        pH = pH + pH_step
        pH_values_ASOR_wt[i] = dep_funct.pH_dependence_ASOR(pH)
        pH_values_ClC[i] = dep_funct.pH_dependence_ClC(pH)
        pH_values_ASOR_mutant[i] = dep_funct.pH_dependence_ASOR(
            pH, pH_k2=1.0, pH_half=7.4
        )

    for i in range(len(voltage_axis)):
        U_values_ASOR[i] = dep_funct.v_dependence_ASOR(voltage_axis[i])
        U_values_ClC[i] = dep_funct.v_dependence_ClC(voltage_axis[i])

    titles = ["ClC", "ASOR", "", ""]
    fig = make_subplots(rows=2, cols=2, subplot_titles=titles)

    position = {"row": 1, "col": 1}
    fig.add_trace(
        go.Scatter(
            x=pH_axis,
            y=pH_values_ClC,
            name=titles.pop(0),
        ),
        **position,
    )
    fig.update_xaxes(title_text="pH", **position)
    fig.update_yaxes(title_text="pH-dependency function", **position)

    position = {"row": 1, "col": 2}
    fig.add_trace(
        go.Scatter(x=pH_axis, y=pH_values_ASOR_wt, name="ASOR wild-type"),
        **position,
    )
    fig.add_trace(
        go.Scatter(
            x=pH_axis,
            y=pH_values_ASOR_mutant,
            name="pH-shifted ASOR mutant",
        ),
        **position,
    )
    fig.update_xaxes(title_text="pH", **position)

    position = {"row": 2, "col": 1}
    fig.add_trace(
        go.Scatter(x=voltage_axis * 1000, y=U_values_ClC, name="U values ClC"),
        **position,
    )
    fig.update_yaxes(title_text="U-dependency function", **position)
    fig.update_xaxes(title_text="U, mV", **position)

    position = {"row": 2, "col": 2}
    fig.add_trace(
        go.Scatter(x=voltage_axis * 1000, y=U_values_ASOR, name="U values ASOR"),
        **position,
    )
    fig.update_xaxes(title_text="U, mV", **position)

    fig.update_layout(height=FIG_HEIGHT, width=FIG_WIDTH)

    return fig
