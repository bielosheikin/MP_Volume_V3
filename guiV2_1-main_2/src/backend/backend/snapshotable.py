from abc import ABC


class Snapshotable(ABC):
    FIELDS_TO_SNAPSHOT = ()

    def __init__(self):
        pass

    def make_snapshot(self):
        snapshot = dict()
        for field in self.FIELDS_TO_SNAPSHOT:
            field_val = getattr(self, field)
            if isinstance(field_val, Snapshotable):
                snapshot[field] = field_val.make_snapshot()
            else:
                snapshot[field] = field_val

        return snapshot
