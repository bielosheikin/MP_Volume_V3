from .snapshotable import Snapshotable

class ExteriorConfig:

    DEFAULT_PH = 7.2

    def __init__(self,
                 *,
                 pH: float = None):
        
        self.pH = pH if pH is not None else self.DEFAULT_PH


class Exterior(Snapshotable):
    FIELDS_TO_SNAPSHOT = ('pH',)

    def __init__(self,
                 *,
                 config: ExteriorConfig = None,
                 init_buffer_capacity: float = None,
                 **kwargs):
        super(Exterior, self).__init__(**kwargs)

        self.config = config if config is not None else ExteriorConfig()
        self.pH = self.config.pH
