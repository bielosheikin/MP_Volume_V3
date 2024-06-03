from .constants import IDEAL_GAS_CONSTANT, FARADAY_CONSTANT

from .snapshotable import Snapshotable

from .exterior import ExteriorConfig, Exterior
from .vesicle import VesicleConfig, Vesicle


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
    FIELDS_TO_SNAPSHOT = ('exterior', 'vesicle')

    def __init__(self,
                 *args,
                 config: SimulationConfig = None,
                 **kwargs):
        super(Simulation, self).__init__(*args, **kwargs)
        self.config = config if config is not None else SimulationConfig()
        self.exterior = Exterior(config=self.config.exterior_config,
                                 init_buffer_capacity=self.config.init_buffer_capacity)
        self.vesicle = Vesicle(simulation_config=self.config.vesicle_config,
                               init_buffer_capacity=self.config.init_buffer_capacity)

    def run(self):
        pass