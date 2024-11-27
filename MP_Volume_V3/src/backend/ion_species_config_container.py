from .ion_species import IonSpecies

class IonSpeciesConfig:
    cl = IonSpecies(
        display_name='cl',
        init_vesicle_conc=159 * 1e-3,
        exterior_conc=20 * 1e-3,
        elementary_charge=-1
    )
    h = IonSpecies(
        display_name='h',
        init_vesicle_conc=7.962143411069938 * 1e-5,
        exterior_conc=12.619146889603859 * 1e-5,
        elementary_charge=1
    )
    na = IonSpecies(
        display_name='na',
        init_vesicle_conc=150 * 1e-3,
        exterior_conc=10 * 1e-3,
        elementary_charge=1
    )
    k = IonSpecies(
        display_name='k',
        init_vesicle_conc=5 * 1e-3,
        exterior_conc=140 * 1e-3,
        elementary_charge=1
    )

    @classmethod
    def update_species(cls, species_name, **kwargs):
        """Update the configuration of a specific ion species."""
        if hasattr(cls, species_name):
            species_config = getattr(cls, species_name)
            for key, value in kwargs.items():
                if hasattr(species_config, key):
                    setattr(species_config, key, value)
                else:
                    raise AttributeError(f"{key} is not a valid attribute for {species_name}.")
        else:
            raise AttributeError(f"No ion species named {species_name} exists.")

    @classmethod
    def get_species(cls, species_name):
        """Retrieve a specific ion species configuration."""
        if hasattr(cls, species_name):
            return getattr(cls, species_name)
        else:
            raise AttributeError(f"No ion species named {species_name} exists.")