from .snapshotable import Snapshotable

class ExteriorConfig:
    DEFAULT_CL_CONCENTRATION = 2e-2
    DEFAULT_NA_CONCENTRATION = 1e-2
    DEFAULT_K_CONCENTRATION = 1.4e-1

    DEFAULT_PH = 7.2

    def __init__(self,
                 *,
                 cl_concentration: float = None,
                 na_concentration: float = None,
                 k_concentration: float = None,
                 pH: float = None):
        self.cl_concentration = cl_concentration if cl_concentration is not None else self.DEFAULT_CL_CONCENTRATION
        self.na_concentration = na_concentration if na_concentration is not None else self.DEFAULT_NA_CONCENTRATION
        self.k_concentration = k_concentration if k_concentration is not None else self.DEFAULT_K_CONCENTRATION

        self.pH = pH if pH is not None else self.DEFAULT_PH


class Exterior(Snapshotable):
    FIELDS_TO_SNAPSHOT = (
        'cl_concentration',
        'na_concentration',
        'k_concentration',
        'pH',
        'free_h_concentration',
        'total_h_concentration'
    )

    def __init__(self,
                 *args,
                 config: ExteriorConfig = None,
                 init_buffer_capacity: float = None,
                 **kwargs):
        super(Exterior, self).__init__(*args, **kwargs)

        self.config = config if ExteriorConfig is not None else ExteriorConfig()
        self.cl_concentration = self.config.cl_concentration
        self.na_concentration = self.config.na_concentration
        self.k_concentration = self.config.k_concentration

        self.pH = self.config.pH
        self.free_h_concentration = 10**(-self.pH)

        # To be filled from outside by Simulation
        self.total_h_concentration = self.free_h_concentration / init_buffer_capacity
