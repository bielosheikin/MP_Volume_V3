from .trackable import Trackable
from .ion_channels import IonChannel
from .flux_calculation_parameters import FluxCalculationParameters


class IonSpecies(Trackable):
    TRACKABLE_FIELDS = ('vesicle_conc', 'vesicle_amount')

    def __init__(self,
                 *,
                 init_vesicle_conc: float,
                 exterior_conc: float,
                 elementary_charge: float,
                 **kwargs
                 ):
        super(IonSpecies, self).__init__(**kwargs)
        self.init_vesicle_conc = init_vesicle_conc
        self.exterior_conc = exterior_conc
        self.elementary_charge = elementary_charge
        self.channels = []
        self.vesicle_conc = self.init_vesicle_conc
        self.vesicle_amount = None
    
    def add_channel(self,
                    *,
                    channel: IonChannel = None):
        self.channels.append(channel)

    def remove_channel(self,
                       *,
                       channel: IonChannel = None):
        if channel in self.channels:
            self.channels.remove(channel)


    def compute_total_flux(self,
                           *,
                           flux_calculation_parameters: FluxCalculationParameters = None
                           ):

        total_flux = 0.0
        for channel in self.channels:
            flux = channel.compute_flux(exterior_conc=self.exterior_conc,
                                        vesicle_conc=self.vesicle_conc,
                                        flux_calculation_parameters=flux_calculation_parameters)
            total_flux += flux
        
        return total_flux