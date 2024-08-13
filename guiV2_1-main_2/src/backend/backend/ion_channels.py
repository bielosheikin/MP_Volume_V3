from math import exp, log

from abc import ABC

from typing import Type

from .snapshotable import Snapshotable
from .vesicle import Vesicle
from .exterior import Exterior


class IonChannelConfig(Snapshotable):
    FIELDS_TO_SNAPSHOT = ('conductance',)

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

    CONFIG_CLASS = None


    def __init__(self,
                 *,
                 config: IonChannelConfig = None,
                 nernst_constant: float = None,
                 ):
        super(IonChannel, self).__init__()
        self.config = config if config is not None else self.CONFIG_CLASS()

        self.conductance = self.config.conductance

        self.nernst_constant = nernst_constant

        self.flux = None
        self.nernst_potential = None

    def get_free_hydrogen_concentrations(self,
                                        vesicle: Vesicle,
                                        exterior: Exterior
                                        ):

        """Retrieve the free hydrogen concentrations from vesicle and exterior."""
        return vesicle.free_h_conc, exterior.free_h_concentration


    def compute_log_term(self, 
                         exterior_conc: float, 
                         vesicle_conc: float,
                         vesicle: Vesicle = None, 
                         exterior: Exterior = None,
                         include_hydrogen: bool = False,
                         exterior_exp: float = 1.0,
                         vesicle_exp: float = 1.0
                         ) -> float:
        
        """Compute the logarithmic term in the Nernst equation, with optional hydrogen ion factor."""
        log_term = (exterior_conc**exterior_exp) / (vesicle_conc**vesicle_exp)

        if include_hydrogen and vesicle and exterior:
            vesicle_h, exterior_h = self.get_free_hydrogen_concentrations(vesicle, exterior)
            log_term *= vesicle_h / exterior_h

        return log(log_term)

    def compute_nernst_potential(self,
                                 voltage: float = None,
                                 exterior_conc: float = None,
                                 vesicle_conc: float = None,
                                 vesicle: Vesicle = None,
                                 exterior: Exterior = None
                                 ):
        
        """Compute the Nernst potential with a customizable log term."""
        log_term = self.compute_log_term(exterior_conc,
                                         vesicle_conc,
                                         vesicle=vesicle,
                                         exterior=exterior)
        
        self.nernst_potential = (self.VOLTAGE_MULTIPLIER * voltage
               + self.NERNST_MULTIPLIER * self.nernst_constant * log_term)

        return self.nernst_potential - self.VOLTAGE_SHIFT

    def compute_flux(self,
                     voltage: float = None,
                     exterior_conc: float = None,
                     vesicle_conc: float = None,
                     area: float = None,
                     vesicle: Vesicle = None,
                     exterior: Exterior = None
                     ):
        self.nernst_potential = self.compute_nernst_potential(voltage=voltage,
                                                              exterior_conc=exterior_conc,
                                                              vesicle_conc=vesicle_conc,
                                                              vesicle=vesicle,
                                                              exterior=exterior)

        return self.FLUX_MULTIPLIER * self.nernst_potential * self.conductance * area


class DependantIonChannel(IonChannel):
    FIELDS_TO_SNAPSHOT = IonChannel.FIELDS_TO_SNAPSHOT + ('pH_dependence', 'voltage_dependence', 'time_dependence')

    def __init__(self,
                 **kwargs):
        super(DependantIonChannel, self).__init__(**kwargs)
        self.pH_exponent = None
        self.half_act_pH = None
        self.voltage_exponent = None
        self.half_act_voltage = None
        self.time_exponent = None
        self.half_act_time = None
        self.time_dependence = None
        self.pH_dependence = None
        self.voltage_dependence = None

    def compute_pH_dependence(self,
                              pH: float = None):
        if self.pH_exponent is None or self.half_act_pH is None:
            raise ValueError("pH dependence parameters are not set.")

        self.pH_dependence = 1.0 / (1.0 + exp(self.pH_exponent * (pH - self.half_act_pH)))
        return self.pH_dependence

    def compute_voltage_dependence(self,
                                   voltage: float = None):
        if self.voltage_exponent is None or self.half_act_voltage is None:
            raise ValueError("Voltage dependence parameters are not set.")
        self.voltage_dependence = 1.0 / (1.0 + exp(self.voltage_exponent * (voltage - self.half_act_voltage)))
        return self.voltage_dependence
    
    def compute_time_dependence(self,
                                t: float = None
                                ):
        if self.time_exponent is None or self.half_act_time is None:
            raise ValueError("Time dependence parameters are not set")
        self.time_dependence = 1.0 / (1.0 + exp(self.time_exponent * (self.half_act_time - t)))
        return self.time_dependence
    
    def compute_flux(self,
                     voltage: float = None,
                     exterior_conc: float = None,
                     vesicle_conc: float = None,
                     area: float = None,
                     vesicle: Vesicle = None,
                     exterior: Exterior = None,
                     t: float = None
                     ):
        

        base_flux = IonChannel.compute_flux(self, 
                                            voltage=voltage,
                                            exterior_conc=exterior_conc,
                                            vesicle_conc=vesicle_conc,
                                            area=area,
                                            vesicle=vesicle,
                                            exterior=exterior)
        flux = base_flux

        if self.pH_dependence is not None:
            flux *= self.pH_dependence

        if self.voltage_dependence is not None:
            flux *=self.voltage_dependence

        if t is not None and self.time_dependence is not None:
            self.compute_time_dependence(t)
            flux *= self.time_dependence
        
        return flux
    

class ASORChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 8e-5

    ALLOWED_CHANNEL_TYPES = ('wt', 'mt', 'none')
    DEFAULT_CHANNEL_TYPE = 'wt'

    ALLOWED_VOLTAGE_DEPS = ('yes', 'no')
    DEFAULT_VOLTAGE_DEP = 'yes'


    def __init__(self,
                 *,
                 channel_type: str = None,
                 voltage_dep: str = None,
                 **kwargs):
        super(ASORChannelConfig, self).__init__(**kwargs)

        channel_type = channel_type if channel_type is not None else self.DEFAULT_CHANNEL_TYPE
        voltage_dep = voltage_dep if voltage_dep is not None else self.DEFAULT_VOLTAGE_DEP

        assert channel_type in self.ALLOWED_CHANNEL_TYPES
        assert voltage_dep in self.ALLOWED_VOLTAGE_DEPS

        self.channel_type = channel_type
        self.voltage_dep = voltage_dep


class ASORChannel(DependantIonChannel):
    
    CONFIG_CLASS = ASORChannelConfig

    def __init__(self,
                 **kwargs):
        super(ASORChannel, self).__init__(**kwargs)

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

        match self.config.voltage_dep:
            case 'yes':
                self.voltage_exponent = 80.0
                self.half_act_voltage = -4e-2
            case 'no':
                self.pH_exponent = 0.0
                self.half_act_pH = 0.0
            case _:
                raise ValueError(f'Wrong voltage dependence: {self.config.voltage_dep}')
    

class TPCChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 2e-6

    def __init__(self,
                 **kwargs):
        super(TPCChannelConfig, self).__init__(**kwargs)

class TPCChannel(IonChannel):

    CONFIG_CLASS = TPCChannelConfig

    def __init__(self,
                 **kwargs):
        super(TPCChannel, self).__init__(**kwargs)


class KChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 0.0

    def __init__(self,
                 **kwargs):
        super(KChannelConfig, self).__init__(**kwargs)

class KChannel(IonChannel):

    CONFIG_CLASS = KChannelConfig

    def __init__(self,
                 **kwargs):
        super(KChannel, self).__init__(**kwargs)


class CLCChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 1e-7

    ALLOWED_DEP_TYPES = ('yes', 'no')
    DEFAULT_DEP_TYPE = 'yes'

    def __init__(self,
                 *,
                 dep_type: str = None,
                 **kwargs):
        super(CLCChannelConfig, self).__init__(**kwargs)
        
        self.dep_type = dep_type if dep_type is not None else self.DEFAULT_DEP_TYPE
        assert self.dep_type in self.ALLOWED_DEP_TYPES


class CLCChannel(DependantIonChannel):
    
    CONFIG_CLASS = CLCChannelConfig

    def __init__(self,
                 **kwargs):
        super(CLCChannel, self).__init__(**kwargs)

        match self.config.dep_type:
            case 'yes':
                self.pH_exponent = -1.5
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
    
    def compute_log_term(self,
                         exterior_conc: float,
                         vesicle_conc: float,
                         vesicle: Vesicle = None,
                         exterior: Exterior = None,
                         include_hydrogen: bool = True
                         ) -> float:

        return super().compute_log_term(exterior_conc, 
                                        vesicle_conc,
                                        vesicle=vesicle, exterior=exterior,
                                        include_hydrogen=include_hydrogen,
                                        exterior_exp=2.0,
                                        vesicle_exp=2.0
                                        )


class NHEChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 0.0

    def __init__(self,
                 **kwargs):
        super(NHEChannelConfig, self).__init__(**kwargs)


class NHEChannel(IonChannel):

    CONFIG_CLASS = NHEChannelConfig

    def __init__(self,
                 **kwargs):
        super(NHEChannel, self).__init__(**kwargs)

    def compute_log_term(self,
                         exterior_conc: float,
                         vesicle_conc: float,
                         vesicle: Vesicle = None,
                         exterior: Exterior = None,
                         include_hydrogen: bool = True
                         ) -> float:
        
        return super().compute_log_term(exterior_conc, 
                                        vesicle_conc,
                                        vesicle=vesicle, 
                                        exterior=exterior,
                                        include_hydrogen=include_hydrogen
                                        )


class VATPaseChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 8e-9

    def __init__(self,
                 **kwargs):
        super(VATPaseChannelConfig, self).__init__(**kwargs)

class VATPaseChannel(DependantIonChannel):

    CONFIG_CLASS = VATPaseChannelConfig

    def __init__(self,
                 **kwargs):
        super(VATPaseChannel, self).__init__(**kwargs)

        self.time_exponent = 0
        self.half_act_time = 0


class HLeakChannelConfig(IonChannelConfig):
    DEFAULT_CONDUCTANCE = 1.6e-8

    def __init__(self,
                 **kwargs):
        super(HLeakChannelConfig, self).__init__(**kwargs)

class HLeakChannel(IonChannel):
    
    CONFIG_CLASS = HLeakChannelConfig

    def __init__(self,
                 **kwargs):
        super(HLeakChannel, self).__init__(**kwargs)

    def compute_log_term(self,
                         exterior_conc, 
                         vesicle_conc, 
                         vesicle: Vesicle = None, 
                         exterior: Exterior = None,
                         include_hydrogen: bool = True
                         ) -> float:
        
        if include_hydrogen and vesicle and exterior:
            vesicle_h, exterior_h = self.get_free_hydrogen_concentrations(vesicle, exterior)
            log_term = exterior_h / vesicle_h
        else:
            raise ValueError("Hydrogen concentration data must be provided for HLeakChannel")
        
        return log(log_term)    