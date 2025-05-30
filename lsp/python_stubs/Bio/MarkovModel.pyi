from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning, MissingPythonDependencyError as MissingPythonDependencyError
from _typeshed import Incomplete

logaddexp: Incomplete

def itemindex(values): ...

VERY_SMALL_NUMBER: float
LOG0: Incomplete

class MarkovModel:
    states: Incomplete
    alphabet: Incomplete
    p_initial: Incomplete
    p_transition: Incomplete
    p_emission: Incomplete
    def __init__(self, states, alphabet, p_initial: Incomplete | None = None, p_transition: Incomplete | None = None, p_emission: Incomplete | None = None) -> None: ...

def load(handle): ...
def save(mm, handle) -> None: ...
def train_bw(states, alphabet, training_data, pseudo_initial: Incomplete | None = None, pseudo_transition: Incomplete | None = None, pseudo_emission: Incomplete | None = None, update_fn: Incomplete | None = None): ...

MAX_ITERATIONS: int

def train_visible(states, alphabet, training_data, pseudo_initial: Incomplete | None = None, pseudo_transition: Incomplete | None = None, pseudo_emission: Incomplete | None = None): ...
def find_states(markov_model, output): ...
