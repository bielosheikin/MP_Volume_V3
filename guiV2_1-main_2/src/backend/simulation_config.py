from copy import deepcopy
from typing import Self

import numpy as np

# --- Constants ---

# mapping of chemicals to numerical indices used to index the state vectors
Cl_idx = 0
Na_idx = 1
H_idx = 2
K_idx = 3

number_of_ions = 4
# T = 1000 # total time (in seconds)
# dt = 0.001 # duration of integration window (in seconds). Has to be 0.001 or less to make integration stable

RT = 2578.5871  # gas constant x temperature in Kelvin
F = 96485.0  # Faraday constant
c_spec = 0.01  # membrane capacitance per surface area F/(m**2)
buffer_capacity_t0 = 5.0 * 1e-4  # initial buffer capacity
U0 = 40 * 1e-3  # membrane potential in V

r = 1.3e-6  # radius of the vesicle in m
V0 = (4.0 / 3.0) * np.pi * (r**3)  # Initial volume of the vesicle in m3
A0 = 4.0 * np.pi * (r**2)  # initial surface area of the vesicle in m2
A_from_V_const = (36.0 * np.pi) ** (
    1 / 3
)  # constant necessary for calculation surface area from the volume

C0 = A0 * c_spec  # inititial membrane capacitance

Cl_o_concentration = 20 * 1e-3  # external Cl concentraion in M

Na_o_concentration = 10 * 1e-3  # external Na concentraion in M
Na_i_concentration = 150 * 1e-3  # internal Na concentraion in M

K_i_concentration = 5 * 1e-3  # internal K concentraion in M
K_o_concentration = 140 * 1e-3  # external Cl concentraion in M

pH_o = 7.2  # external pH
pH_i = 7.4  # internal pH

hfree_o_concentration = 10 ** (-pH_o)  # concentration of free external protons in M
hfree_i_concentration = 10 ** (-pH_i)  # concentration of free internal protons in M

htotal_o_concentration = (
    hfree_o_concentration / buffer_capacity_t0
)  # concentration of total external protons in M
htotal_i_concentration = (
    hfree_i_concentration / buffer_capacity_t0
)  # concentration of total intenral protons in M

htotal_i_amount = (
    htotal_i_concentration * V0 * 1000
)  # amount of total internal protons in moles


class SimulationConfig:
    def __init__(self) -> None:
        self.parameters = dict()
        self.additional_parameters = dict()

    def update(self, input_fields: dict) -> Self:
        params = deepcopy(input_fields)

        # Handle ASOR params
        ASOR_args = {}
        match params["ASOR"]["pH_dep"]:
            case "wt":
                ASOR_args["ASOR_pH_k2"] = 3.0
                ASOR_args["ASOR_pH_half"] = 5.4
            case "mt":
                ASOR_args["ASOR_pH_k2"] = 1.0
                ASOR_args["ASOR_pH_half"] = 7.4
            case _:
                ASOR_args["ASOR_pH_k2"] = 0.0
                ASOR_args["ASOR_pH_half"] = 0.0

        match params["ASOR"]["U_dep"]:
            case True:
                ASOR_args["ASOR_U_k2"] = 80.0
                ASOR_args["ASOR_U_half"] = -40 * 1e-3
            case False:
                ASOR_args["ASOR_U_k2"] = 0.0
                ASOR_args["ASOR_U_half"] = 0.0
            case _:
                raise ValueError("Wrong ASOR U_dep param value. Could be True or False")
        # Handle CLC params
        CLC_args = {}
        match params["Other"]["CLC_dep"]:
            case True:
                CLC_args["CLC_pH_k2"] = 1.5
                CLC_args["CLC_pH_half"] = 5.5
                CLC_args["CLC_U_k2"] = 80.0
                CLC_args["CLC_U_half"] = -40 * 1e-3
            case False:
                CLC_args["CLC_pH_k2"] = 0
                CLC_args["CLC_pH_half"] = 0
                CLC_args["CLC_U_k2"] = 0
                CLC_args["CLC_U_half"] = 0
            case _:
                raise ValueError("Wrong CLC_dep  param value. Could be True or False")

        # Handle Cl_i param

        match params["Other"]["Cl_i"]:
            case "high":
                Cl_i_concentration = 159 * 1e-3
            case "low":
                Cl_i_concentration = 9 * 1e-3
            case "zero":
                Cl_i_concentration = 1 * 1e-3
            case _:
                raise ValueError(
                    "Wrong CL_i param value.Could be `high`, `low` or `zero`"
                )
        (
            X_amount,
            external_ions_concentrations,
            internal_ions_amounts,
            internal_ions_concentrations,
            Sum_initial_amounts,
        ) = self.initialize_internal_concentrations(Cl_i_concentration)

        self.parameters = {
            "dt": params["Time"]["dt"],
            "T": params["Time"]["T"],
            "G": params["G"],
            "external_ions_concentrations": external_ions_concentrations,
            "A_from_V_const": A_from_V_const,
            "X_amount": X_amount,
            "buffer_capacity_t0": buffer_capacity_t0,
            "V_t0": V0,
            "c_spec": c_spec,
            "RT": RT,
            "F": F,
            "pH_i": pH_i,
            "U0": U0,
            "A0": A0,
            "C0": C0,
            "Sum_initial_amounts": Sum_initial_amounts,
            "ASOR_pH_k2": ASOR_args["ASOR_pH_k2"],
            "ASOR_pH_half": ASOR_args["ASOR_pH_half"],
            "ASOR_U_k2": ASOR_args["ASOR_U_k2"],
            "ASOR_U_half": ASOR_args["ASOR_U_half"],
            "CLC_pH_k2": CLC_args["CLC_pH_k2"],
            "CLC_pH_half": CLC_args["CLC_pH_half"],
            "CLC_U_k2": CLC_args["CLC_U_k2"],
            "CLC_U_half": CLC_args["CLC_U_half"],
        }

        self.additional_parameters = {
            "internal_ions_amounts": internal_ions_amounts,
            "internal_ions_concentrations": internal_ions_concentrations,
        }

        return self

    def initialize_internal_concentrations(self, Cl_i_concentration=159 * 1e-3):
        # Cl_i_concentration= 159*1e-3 # internal Cl concentraion in M
        # Cl_i_concentration=1*1e-3 # absent internal Cl concentration in M
        # Cl_i_concentration= 9*1e-3 # Cl replacement condition (internal Cl concentration in M)

        Q0 = U0 * C0  # initial total charge
        X_amount = (Q0 / F) - (
            (
                Na_i_concentration
                + K_i_concentration
                + htotal_i_concentration
                - Cl_i_concentration
            )
            * V0
            * 1000
        )  # initial amount of unaccouted ions in moles
        X_concentration = X_amount / (
            V0 * 1000
        )  # initial concentration of unaccounted ions in moles

        internal_ions_amounts = [
            Cl_i_concentration * V0 * 1000,
            Na_i_concentration * V0 * 1000,
            htotal_i_concentration * V0 * 1000,
            K_i_concentration * V0 * 1000,
        ]  # vector of amounts of ions in moles
        external_ions_concentrations = [
            Cl_o_concentration,
            Na_o_concentration,
            htotal_o_concentration,
            K_o_concentration,
        ]  # vector of concentrations of external ions
        internal_ions_concentrations = [
            Cl_i_concentration,
            Na_i_concentration,
            htotal_i_concentration,
            K_i_concentration,
        ]  # vector of concentrations of internal ions

        Sum_initial_amounts = (
            internal_ions_amounts[Cl_idx]
            + internal_ions_amounts[Na_idx]
            + abs(X_amount)
            + internal_ions_amounts[K_idx]
        )  # sum of amounts of all ions

        return (
            X_amount,
            external_ions_concentrations,
            internal_ions_amounts,
            internal_ions_concentrations,
            Sum_initial_amounts,
        )
