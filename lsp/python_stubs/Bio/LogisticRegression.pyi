from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning, MissingPythonDependencyError as MissingPythonDependencyError
from _typeshed import Incomplete

class LogisticRegression:
    beta: Incomplete
    def __init__(self) -> None: ...

def train(xs, ys, update_fn: Incomplete | None = None, typecode: Incomplete | None = None): ...
def calculate(lr, x): ...
def classify(lr, x): ...
