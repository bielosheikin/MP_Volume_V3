from .constants import IDEAL_GAS_CONSTANT, FARADAY_CONSTANT, VOLUME_TO_AREA_CONSTANT

from .trackable import Trackable

from .exterior import ExteriorConfig, Exterior
from .vesicle import VesicleConfig, Vesicle
from .ion_species import IonSpecies
from .ion_channels import IonChannel
from .flux_calculation_parameters import FluxCalculationParameters
from .histories_storage import HistoriesStorage
from math import log10


class SimulationConfig:
    DEFAULT_TIME_STEP = 0.001
    DEFAULT_TOTAL_TIME = 100.0

    DEFAULT_TEMPERATURE = 2578.5871 / IDEAL_GAS_CONSTANT
    DEFAULT_INIT_BUFFER_CAPACITY = 5e-4

    def __init__(self,
                 *,
                 time_step: float = None,
                 total_time: float = None,
                 temperature: float = None,
                 init_buffer_capacity: float = None):
        
        self.time_step = time_step if time_step is not None else self.DEFAULT_TIME_STEP
        self.total_time = total_time if total_time is not None else self.DEFAULT_TOTAL_TIME
        self.temperature = temperature if temperature is not None else self.DEFAULT_TEMPERATURE
        self.init_buffer_capacity = init_buffer_capacity if init_buffer_capacity is not None else self.DEFAULT_INIT_BUFFER_CAPACITY
        

class Simulation(Trackable):
    TRACKABLE_FIELDS = ('buffer_capacity','time')

    def __init__(self,
                 *,
                 config: SimulationConfig = None,
                 display_name: str = 'simulation',
                 **kwargs):
        super(Simulation, self).__init__(display_name=display_name, **kwargs)

        self.config = config if config is not None else SimulationConfig()
        self.time = 0.0
        self.exterior = None
        self.vesicle = None
        self.hydrogen_species = None
        self.buffer_capacity = self.config.init_buffer_capacity
        self.all_species = []  
        self.histories = HistoriesStorage()
        self.nernst_constant = self.config.DEFAULT_TEMPERATURE * IDEAL_GAS_CONSTANT / FARADAY_CONSTANT

        self.histories.register_object(self)

    def add_vesicle(self,
                    vesicle_obj: Vesicle
                    ):
        self.vesicle = vesicle_obj
        self.histories.register_object(self.vesicle)

    def add_exterior(self,
                     exterior_obj: Exterior
                     ):
        self.exterior = exterior_obj
        self.histories.register_object(self.exterior)

    def add_hydrogen_species(self, 
                            hydrogen_obj: IonSpecies
                            ):
        self.hydrogen_species = hydrogen_obj

    def add_ion_species(self, 
                        species_obj: IonSpecies
                        ):
        self.all_species.append(species_obj)
        self.histories.register_object(species_obj)
    
    def add_channel(self, 
                    species_obj : IonSpecies, 
                    channel_obj: IonChannel
                    ):
        species_obj.add_channel(channel = channel_obj)
        self.histories.register_object(channel_obj)
    
    def get_Flux_Calculation_Parameters(self):
        flux_calculation_parameters = FluxCalculationParameters()      

        flux_calculation_parameters.voltage = self.vesicle.config.voltage
        flux_calculation_parameters.pH = self.vesicle.pH
        flux_calculation_parameters.area = self.vesicle.area
        flux_calculation_parameters.time = self.time
        flux_calculation_parameters.vesicle_hydrogen_free = self.hydrogen_species.vesicle_conc * self.buffer_capacity
        flux_calculation_parameters.exterior_hydrogen_free = self.hydrogen_species.exterior_conc * self.config.init_buffer_capacity
        flux_calculation_parameters.nernst_constant = self.nernst_constant

        return flux_calculation_parameters
    
    def get_unaccounted_ion_amount(self):

        return ((self.vesicle.init_charge / FARADAY_CONSTANT) - 
                (sum(ion.elementary_charge * ion.init_vesicle_conc for ion in self.all_species)) * 1000 * self.vesicle.init_volume)
    
    def update_volume(self):
        self.vesicle.volume = (self.vesicle.init_volume * 
                               (sum(ion.vesicle_conc for ion in self.all_species if ion.display_name != 'h') + 
                                self.get_unaccounted_ion_amount()) /
                                (sum(ion.init_vesicle_conc for ion in self.all_species if ion.display_name != 'h') +
                                self.get_unaccounted_ion_amount())
                              )

    def update_area(self):
        self.vesicle.area = (VOLUME_TO_AREA_CONSTANT * self.vesicle.volume**(2/3))

    def update_charge(self):
        self.vesicle.charge = ((sum(ion.elementary_charge * ion.vesicle_amount for ion in self.all_species) +
                               self.get_unaccounted_ion_amount()) *
                               FARADAY_CONSTANT
                              )
    
    def update_capacitance(self):
        self.vesicle.capacitance = self.vesicle.area * self.vesicle.config.specific_capacitance
    
    def update_voltage(self):
        self.vesicle.voltage = self.vesicle.charge / self.vesicle.capacitance

    def update_buffer(self):
        self.buffer_capacity = self.config.init_buffer_capacity * self.vesicle.volume / self.vesicle.init_volume

    def update_pH(self):
        self.vesicle.pH = -log10(self.hydrogen_species.vesicle_conc * self.buffer_capacity)

    def set_ion_amounts(self):
        for ion in self.all_species:
            ion.vesicle_amount = ion.vesicle_conc * 1000 * self.vesicle.volume

    def update_ion_amounts(self, fluxes):
        for ion, flux in zip(self.all_species, fluxes):
            ion.vesicle_amount += flux * self.config.time_step
            
            if ion.vesicle_amount < 0:
                ion.vesicle_amount = 0
                print(f"Warning: {ion.display_name} ion amount fell below zero and has been reset to zero.")


    def update_vesicle_concentrations(self):
        for ion in self.all_species:
            ion.vesicle_conc = ion.vesicle_amount / (1000 * self.vesicle.volume)

    def update_simulation_state(self):
        self.update_volume()
        self.update_charge()
        self.update_area()
        self.update_capacitance()
        self.update_voltage()
        self.update_buffer()
        self.update_pH()

    
    def run_one_iteration(self): 

        self.update_simulation_state()
        self.update_vesicle_concentrations()

        flux_calculation_parameters = self.get_Flux_Calculation_Parameters()
        fluxes = [ion.compute_total_flux(flux_calculation_parameters=flux_calculation_parameters) for ion in self.all_species]

        self.histories.update_histories()
        
        self.update_ion_amounts(fluxes)

        self.time += self.config.time_step

    def run(self):

        if self.hydrogen_species is None:
            raise ValueError("Hydrogen species must be provided before running the simulation.")
        
        self.set_ion_amounts()

        while self.time < self.config.total_time:
            self.run_one_iteration()
        
        return self.histories
