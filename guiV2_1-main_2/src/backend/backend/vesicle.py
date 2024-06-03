import math

from .snapshotable import Snapshotable

from .constants import FARADAY_CONSTANT


class VesicleConfig:
    DEFAULT_SPECIFIC_CAPACITANCE = 1e-2
    DEFAULT_INIT_VOLTAGE = 4e-2

    DEFAULT_INIT_RADIUS = 1.3e-6

    ALLOWED_INIT_CL_CONC_TYPES = ('high', 'low', 'zero')
    DEFAULT_INIT_CL_CONC_TYPE = 'high'
    DEFAULT_INIT_NA_CONC = 1.5e-1
    DEFAULT_INIT_K_CONC = 5e-3

    DEFAULT_INIT_PH = 7.4

    def __init__(self,
                 *,
                 specific_capacitance: float = None,
                 init_voltage: float = None,
                 init_radius: float = None,
                 init_cl_conc_type: str = None,
                 init_na_conc: float = None,
                 init_k_conc: float = None,
                 init_pH: float = None):
        self.specific_capacitance = specific_capacitance if specific_capacitance is not None else self.DEFAULT_SPECIFIC_CAPACITANCE
        self.init_voltage = init_voltage if init_voltage is not None else self.DEFAULT_INIT_VOLTAGE

        self.init_radius = init_radius if init_radius is not None else self.DEFAULT_INIT_RADIUS

        if init_cl_conc_type is not None:
            assert init_cl_conc_type in self.ALLOWED_INIT_CL_CONC_TYPES
            self.init_cl_conc_type = init_cl_conc_type
        else:
            self.init_cl_conc_type = self.DEFAULT_INIT_CL_CONC_TYPE
        self.init_na_conc = init_na_conc if init_na_conc is not None else self.DEFAULT_INIT_NA_CONC
        self.init_k_conc = init_k_conc if init_k_conc is not None else self.DEFAULT_INIT_K_CONC

        self.init_pH = init_pH if init_pH is not None else self.DEFAULT_INIT_PH


class Vesicle(Snapshotable):
    FIELDS_TO_SNAPSHOT = (
        'cl_conc',
        'na_conc',
        'k_conc',
        'pH',
        'free_h_conc',
        'total_h_conc',
        'unacc_ion_amount',
        'unacc_ion_conc',
    )

    def __init__(self,
                 *args,
                 config: VesicleConfig = None,
                 init_buffer_capacity: float = None,
                 **kwargs):
        super(Vesicle, self).__init__(*args, **kwargs)

        self.config = config if config is not None else VesicleConfig()
        self.init_volume = (4 / 3) * math.pi * (self.config.init_radius ** 3)
        self.init_area = 4.0 * math.pi * (self.config.init_radius ** 2)

        self.init_capacitance = self.init_area * self.config.specific_capacitance

        self.init_charge = self.config.init_voltage * self.init_capacitance

        self.cl_conc = None
        match self.config.init_cl_conc_type:
            case 'high':
                self.cl_conc = 1.59e-1
            case 'low':
                self.cl_conc = 9e-3
            case 'zero':
                self.cl_conc = 1e-3
            case _:
                raise ValueError(f'Wrong value of init_cl_conc_type: {self.config.init_cl_conc_type}')
        self.na_conc = self.config.init_na_conc
        self.k_conc = self.config.init_k_conc

        self.pH = self.config.init_pH
        self.free_h_conc = 10 ** (-self.pH)

        # To be filled from outside by Simulation
        self.total_h_conc = self.free_h_conc / init_buffer_capacity

        self.unacc_ion_amount = (self.init_charge / FARADAY_CONSTANT) - 1000 * self.init_volume * (
                self.na_conc + self.k_conc + self.total_h_conc - self.cl_conc
        )
        self.unacc_ion_conc = self.unacc_ion_amount / (1000 * self.init_volume)
