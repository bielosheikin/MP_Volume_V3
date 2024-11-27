from .ion_channels import IonChannelConfig

class IonChannelConfig:

    asor = IonChannelConfig(
        conductance=8e-5,
        channel_type='wt',
        voltage_dep='yes',
        dependence_type='voltage_and_pH',
        voltage_multiplier=1,
        nernst_multiplier=1,
        voltage_shift=0,
        flux_multiplier=1,
        allowed_primary_ion='cl',
        primary_exponent=1,
        display_name='ASOR_Config'
    )
    clc = IonChannelConfig(
        conductance=1e-7,
        channel_type='clc',
        voltage_dep='yes',
        dependence_type='voltage_and_pH',
        voltage_multiplier=1,
        nernst_multiplier=1 / 3,
        voltage_shift=0,
        flux_multiplier=2,
        allowed_primary_ion='cl',
        allowed_secondary_ion='h',
        primary_exponent=2,
        secondary_exponent=1,
        use_free_hydrogen=True,
        display_name='CLC_Config'
    )
    tpc = IonChannelConfig(
        conductance=2e-6,
        dependence_type=None,
        voltage_multiplier=-1,
        nernst_multiplier=1,
        voltage_shift=0,
        flux_multiplier=1,
        allowed_primary_ion='na',
        primary_exponent=1,
        display_name='TPC_Config'
    )
    nhe = IonChannelConfig(
        conductance=0.0,
        dependence_type=None,
        voltage_multiplier=0,
        nernst_multiplier=1,
        voltage_shift=0,
        flux_multiplier=1,
        allowed_primary_ion='na',
        allowed_secondary_ion='h',
        custom_nernst_constant=1,
        primary_exponent=1,
        secondary_exponent=1,
        use_free_hydrogen=True,
        display_name='NHE_Config'
    )
    vatpase = IonChannelConfig(
        conductance=8e-9,
        dependence_type='time',
        voltage_multiplier=1,
        nernst_multiplier=-1,
        voltage_shift=0.27,
        flux_multiplier=-1,
        allowed_primary_ion='h',
        primary_exponent=1,
        display_name='VATPase_Config'
    )
    clc_h = IonChannelConfig(
        conductance=1e-7,
        channel_type='clc',
        voltage_dep='yes',
        dependence_type='voltage_and_pH',
        voltage_multiplier=1,
        nernst_multiplier=1 / 3,
        voltage_shift=0,
        flux_multiplier=-1,
        allowed_primary_ion='cl',
        allowed_secondary_ion='h',
        primary_exponent=2,
        secondary_exponent=1,
        use_free_hydrogen=True,
        display_name='CLC_Config_H'
    )
    nhe_h = IonChannelConfig(
        conductance=0.0,
        dependence_type=None,
        voltage_multiplier=0,
        nernst_multiplier=1,
        voltage_shift=0,
        flux_multiplier=-1,
        allowed_primary_ion='na',
        allowed_secondary_ion='h',
        custom_nernst_constant=1,
        primary_exponent=1,
        secondary_exponent=1,
        use_free_hydrogen=True,
        display_name='NHE_Config_H'
    )
    hleak = IonChannelConfig(
        conductance=1.6e-8,
        dependence_type=None,
        voltage_multiplier=-1,
        nernst_multiplier=1,
        voltage_shift=0,
        flux_multiplier=1,
        allowed_primary_ion='h',
        primary_exponent=1,
        use_free_hydrogen=True,
        display_name='HLeak_Config'
    )
    k_channel = IonChannelConfig(
        conductance=0.0,
        dependence_type=None,
        voltage_multiplier=-1,
        nernst_multiplier=1,
        voltage_shift=0,
        flux_multiplier=1,
        allowed_primary_ion='k',
        primary_exponent=1,
        display_name='K_Config'
    )

    @classmethod
    def update_config(cls, channel_name, **kwargs):
        """Update the configuration of a specific channel."""
        if hasattr(cls, channel_name):
            channel_config = getattr(cls, channel_name)
            for key, value in kwargs.items():
                if hasattr(channel_config, key):
                    setattr(channel_config, key, value)
                else:
                    raise AttributeError(f"{key} is not a valid attribute for {channel_name}.")
        else:
            raise AttributeError(f"No channel named {channel_name} exists.")

    @classmethod
    def get_config(cls, channel_name):
        """Retrieve a specific channel configuration."""
        if hasattr(cls, channel_name):
            return getattr(cls, channel_name)
        else:
            raise AttributeError(f"No channel named {channel_name} exists.")