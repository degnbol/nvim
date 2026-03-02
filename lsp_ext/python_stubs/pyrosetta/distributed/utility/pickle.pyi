from _typeshed import Incomplete
from pyrosetta.distributed import with_lock as with_lock

logger: Incomplete

@with_lock
def __cereal_getstate__(self): ...
@with_lock
def __cereal_setstate__(self, state) -> None: ...
def set_cereal_pickleable(klass) -> None: ...
