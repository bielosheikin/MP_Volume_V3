from math import exp, log

from abc import ABC

from typing import Type

from .snapshotable import Snapshotable


class IonChannelConfig:
    DEFAULT_CONDUCTANCE = None

    def __init__(self,
                 *,
                 conductance: float = None):
        self.conductance = conductance if conductance is not None else self.DEFAULT_CONDUCTANCE


class IonChannel(Snapshotable):
    FIELDS_TO_SNAPSHOT = ('flux', 'nernst_potential')
    VOLTAGE_MULTIPLIER = None
    NERNST_MULTIPLIER = None
    VOLTAGE_SHIFT = None

    FLUX_MULTIPLIER = None

    def __init__(self,
                 *args,
                 config_class: Type = None,
                 config: IonChannelConfig = None,
                 nernst_constant: float = None,
                 **kwargs):
        super(IonChannel, self).__init__(*args, **kwargs)
        self.config = config if config is not None else config_class()

        self.conductance = self.config.conductance

        # To be filled by
        self.nernst_constant = nernst_constant

        self.flux = None
        self.nernst_potential = None

    def compute_nernst_potential(self,
                                 voltage: float = None,
                                 exterior_conc: float = None,
                                 vesicle_conc: float = None):
        self.nernst_potential = (self.VOLTAGE_MULTIPLIER * voltage
               + self.NERNST_MULTIPLIER * self.nernst_constant * log(exterior_conc / vesicle_conc))

        return self.nernst_potential - self.VOLTAGE_SHIFT

    def compute_flux(self,
                     voltage: float = None,
                     exterior_conc: float = None,
                     vesicle_conc: float = None,
                     area: float = None):
        self.nernst_potential = self.compute_nernst_potential(voltage=voltage,
                                                              exterior_conc=exterior_conc,
                                                              vesicle_conc=vesicle_conc)

        return self.FLUX_MULTIPLIER * self.nernst_potential * self.conductance * area


class DependantIonChannel(IonChannel):
    FIELDS_TO_SNAPSHOT = IonChannel.FIELDS_TO_SNAPSHOT + ('pH_dependence', 'voltage_dependence')


class ASORChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 8e-5

    ALLOWED_CHANNEL_TYPES = ('wt', 'mt', 'none')
    DEFAULT_CHANNEL_TYPE = 'wt'

    ALLOWED_VOLTAGE_DEPS = ('yes', 'no')
    DEFAULT_VOLTAGE_DEP = 'wt'

    def __init__(self,
                 *args,
                 channel_type: str = None,
                 voltage_dep: str = None,
                 **kwargs):
        super(ASORChannelConfig, self).__init__(*args, **kwargs)
        assert channel_type in self.ALLOWED_CHANNEL_TYPES
        self.channel_type = channel_type if channel_type is not None else self.DEFAULT_CHANNEL_TYPE

        assert voltage_dep in self.ALLOWED_VOLTAGE_DEPS
        self.voltage_dep = voltage_dep if voltage_dep is not None else self.DEFAULT_VOLTAGE_DEP


class ASORChannel(Snapshotable):
    FIELDS_TO_SNAPSHOT = ('pH_dependence', 'voltage_dependence')

    def __init__(self,
                 *args,
                 config: ASORChannelConfig = None,
                 **kwargs):
        super(ASORChannel, self).__init__(*args, **kwargs)

        self.config = config if config is not None else ASORChannelConfig()
        self.pH_exponent = None
        self.half_act_pH = None

        match self.config.channel_type:
            case 'wt':
                self.pH_exponent = 3.0
                self.half_act_pH = 5.4
            case 'mt':
                self.pH_exponent = 1.0
                self.half_act_pH = 7.4
            case 'none':
                self.pH_exponent = 0.0
                self.half_act_pH = 0.0
            case _:
                raise ValueError(f'Wrong channel type: {self.config.channel_type}')

        self.voltage_exponent = None
        self.half_act_voltage = None

        match self.config.voltage_dep:
            case 'yes':
                self.voltage_exponent = 80.0
                self.half_act_voltage = -4e-2
            case 'no':
                self.pH_exponent = 0.0
                self.half_act_pH = 0.0
            case _:
                raise ValueError(f'Wrong voltage dependence: {self.config.voltage_dep}')

        self.pH_dependence = None
        self.voltage_dependence = None

    def compute_pH_dependence(self, pH: float = None):
        self.pH_dependence = 1.0 / (1.0 + exp(self.pH_exponent * (pH - self.half_act_pH)))
        return self.pH_dependence

    def compute_voltage_dependence(self, voltage: float = None):
        self.voltage_dependence = 1.0 / (1.0 + exp(self.voltage_exponent * (voltage - self.half_act_voltage)))

        return self.voltage_dependence


class TPCChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 2e-6

    def __init__(self,
                 *args,
                 **kwargs):
        super(TPCChannelConfig, self).__init__(*args, **kwargs)


class KChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 0.0

    def __init__(self,
                 *args,
                 **kwargs):
        super(KChannelConfig, self).__init__(*args, **kwargs)


class CLCChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 1e-7

    ALLOWED_DEP_TYPES = ('yes', 'no')
    DEFAULT_DEP_TYPE = 'yes'

    def __init__(self,
                 *args,
                 dep_type: str = None,
                 **kwargs):
        super(CLCChannelConfig, self).__init__(*args, **kwargs)
        assert dep_type in self.ALLOWED_DEP_TYPES
        self.dep_type = dep_type if dep_type is not None else self.DEFAULT_DEP_TYPE


class CLCChannel(Snapshotable):
    FIELDS_TO_SNAPSHOT = ('pH_dependece', 'voltage_dependence')

    def __init__(self,
                 *args,
                 config: CLCChannelConfig = None,
                 **kwargs):
        super(CLCChannel, self).__init__(*args, **kwargs)

        self.config = config if config is not None else CLCChannelConfig()
        self.pH_exponent = None
        self.half_act_pH = None
        self.voltage_exponent = None
        self.half_act_voltage = None

        match self.config.dep_type:
            case 'yes':
                self.pH_exponent = 1.5
                self.half_act_pH = 5.5
                self.voltage_exponent = 80.0
                self.half_act_voltage = -4e-2
            case 'no':
                self.pH_exponent = 0.0
                self.half_act_pH = 0.0
                self.voltage_exponent = 0.0
                self.half_act_voltage = 0.0
            case _:
                raise ValueError(f'Wrong dependence type: {self.config.dep_type}')

        self.pH_dependence = None
        self.voltage_dependence = None

    def compute_pH_dependence(self, pH: float = None):
        self.pH_dependence = 1.0 / (1.0 + exp(self.pH_exponent * (pH - self.half_act_pH)))
        return self.pH_dependence

    def compute_voltage_dependence(self, voltage: float = None):
        self.voltage_dependence = 1.0 / (1.0 + exp(self.voltage_exponent * (voltage - self.half_act_voltage)))

        return self.voltage_dependence


class NHEChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 0.0

    def __init__(self,
                 *args,
                 **kwargs):
        super(NHEChannelConfig, self).__init__(*args, **kwargs)


class VATPaseChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 8e-9

    def __init__(self,
                 *args,
                 **kwargs):
        super(VATPaseChannelConfig, self).__init__(*args, **kwargs)


class HLeakChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 1.6e-8

    def __init__(self,
                 *args,
                 **kwargs):
        super(HLeakChannelConfig, self).__init__(*args, **kwargs)
