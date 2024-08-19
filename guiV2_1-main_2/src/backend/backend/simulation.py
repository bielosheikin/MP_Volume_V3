from .constants import IDEAL_GAS_CONSTANT, FARADAY_CONSTANT, VOLUME_TO_AREA_CONSTANT

from .snapshotable import Snapshotable

from .exterior import ExteriorConfig, Exterior
from .vesicle import VesicleConfig, Vesicle
from .ion_species import IonSpecies
from .flux_calculation_parameters import FluxCalculationParameters


class SimulationConfig:
    DEFAULT_TIME_STEP = 0.001
    DEFAULT_TOTAL_TIME = 1000.0

    DEFAULT_TEMPERATURE = 2578.5871 / IDEAL_GAS_CONSTANT
    DEFAULT_INIT_BUFFER_CAPACITY = 5e-4

    def __init__(self,
                 *,
                 time_step: float = None,
                 total_time: float = None,
                 temperature: float = None,
                 exterior_config: ExteriorConfig = None,
                 init_buffer_capacity: float = None,
                 vesicle_config: VesicleConfig = None):
        
        self.time_step = time_step if time_step is not None else self.DEFAULT_TIME_STEP
        self.total_time = total_time if total_time is not None else self.DEFAULT_TOTAL_TIME

        self.temperature = temperature if temperature is not None else self.DEFAULT_TEMPERATURE
        self.exterior_config = exterior_config if exterior_config is not None else ExteriorConfig()
        self.init_buffer_capacity = init_buffer_capacity if init_buffer_capacity is not None else self.DEFAULT_INIT_BUFFER_CAPACITY
        self.vesicle_config = vesicle_config if vesicle_config is not None else VesicleConfig()
        

class Simulation(Snapshotable):
    FIELDS_TO_SNAPSHOT = ('ion_species',)

    def __init__(self,
                 *,
                 config: SimulationConfig = None,
                 **kwargs):
        super(Simulation, self).__init__(**kwargs)

        self.config = config if config is not None else SimulationConfig()
        self.time = 0.0
        self.exterior = Exterior(config=self.config.exterior_config,
                                 init_buffer_capacity=self.config.init_buffer_capacity)
        self.vesicle = Vesicle(config=self.config.vesicle_config,
                               init_buffer_capacity=self.config.init_buffer_capacity)
        self.hydrogene_species = None
        self.buffer_capacity = self.config.init_buffer_capacity
        self.all_species = []  

        self.initialize_hydrogen_species(
            init_vesicle_conc=7.962143411069938e-5,
            exterior_conc=12.619146889603859e-5
        )        

    def initialize_hydrogen_species(self, 
                                    init_vesicle_conc: float, 
                                    exterior_conc: float):
        hydrogen = IonSpecies(
            name='H',
            init_vesicle_conc=init_vesicle_conc,
            exterior_conc=exterior_conc,
            elementary_charge=1
        )
        self.all_species.append(hydrogen)
    
    def get_Flux_Calculation_Parameters(self):
        flux_calculation_parameters = FluxCalculationParameters()

        hydrogen = next((ion for ion in self.all_species if ion.name == 'H'), None)

        if not hydrogen:
            raise ValueError("Hydrogen species not initialized in the simulation.")        

        flux_calculation_parameters.voltage = self.vesicle.config.voltage
        flux_calculation_parameters.pH = self.vesicle.pH
        flux_calculation_parameters.area = self.vesicle.area
        flux_calculation_parameters.time = self.time
        flux_calculation_parameters.vesicle_hydrogene_free = hydrogen.vesicle_conc * self.buffer_capacity
        flux_calculation_parameters.exterior_hydrogene_free = hydrogen.exterior_conc * self.config.init_buffer_capacity
        flux_calculation_parameters.nernst_constant = self.config.DEFAULT_TEMPERATURE * IDEAL_GAS_CONSTANT / FARADAY_CONSTANT

        return flux_calculation_parameters
    
    def get_unaccounted_ion_amount(self):

        return ((self.vesicle.init_charge / FARADAY_CONSTANT) - 
                (sum(ion.elementary_charge * ion.vesicle_conc for ion in self.ion_species) * 1000 * self.vesicle.init_volume))
    
    def update_volume(self):
        self.vesicle.volume = (self.vesicle.init_volume * 
                               (sum(ion.vesicle_conc for ion in self.ion_species if ion.name != 'H') + 
                                self.get_unaccounted_ion_amount()) /
                                (sum(ion.init_vesicle_conc for ion in self.ion_species if ion.name != 'H') +
                                self.get_unaccounted_ion_amount())
                              )

    def update_area(self):
        self.vesicle.area = (VOLUME_TO_AREA_CONSTANT * self.vesicle.volume)

    def update_charge(self):
        self.vesicle.charge = ((sum(ion.elementary_charge * ion.vesicle_conc for ion in self.ion_species if ion.name != 'H') +
                               self.get_unaccounted_ion_amount()) *
                               FARADAY_CONSTANT
                              )
    
    def update_capacitance(self):
        self.vesicle.capacitance = self.vesicle.area * self.vesicle.config.specific_capacitance
    
    def update_voltage(self):
        self.vesicle.config.voltage = self.vesicle.charge / self.vesicle.capacitance

    def update_buffer(self):
        self.buffer_capacity = self.config.init_buffer_capacity * self.vesicle.volume / self.vesicle.init_volume
        
    def add_ion_species(self, 
                        ion_species: IonSpecies = None):
         self.ion_species.append(ion_species)

    def update_ion_concentrations(self, fluxes):
        for ion, flux in zip(self.all_species, fluxes):
   
            ion.vesicle_conc += flux * self.config.time_step

    def run(self):

        if not self.config:
            raise ValueError("No config provided")
        
        result = self.run_simulation(self.config)

        return result
    
    def run_one_iteration(self, config: SimulationConfig):
            
        flux_calculation_parameters = self.get_Flux_Calculation_Parameters()
        fluxes = [ion.compute_total_flux(flux_calculation_parameters=flux_calculation_parameters) for ion in self.all_species]

        self.update_ion_concentrations(fluxes)
        
        # Update everything else: 

        self.time += config.time_step