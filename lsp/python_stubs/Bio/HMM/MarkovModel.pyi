from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning
from Bio.Seq import Seq as Seq
from _typeshed import Incomplete

class MarkovModelBuilder:
    DEFAULT_PSEUDO: int
    initial_prob: Incomplete
    transition_prob: Incomplete
    emission_prob: Incomplete
    transition_pseudo: Incomplete
    emission_pseudo: Incomplete
    def __init__(self, state_alphabet, emission_alphabet) -> None: ...
    def get_markov_model(self): ...
    def set_initial_probabilities(self, initial_prob) -> None: ...
    def set_equal_probabilities(self) -> None: ...
    def set_random_initial_probabilities(self): ...
    def set_random_transition_probabilities(self): ...
    def set_random_emission_probabilities(self): ...
    def set_random_probabilities(self) -> None: ...
    def allow_all_transitions(self) -> None: ...
    def allow_transition(self, from_state, to_state, probability: Incomplete | None = None, pseudocount: Incomplete | None = None) -> None: ...
    def destroy_transition(self, from_state, to_state) -> None: ...
    def set_transition_score(self, from_state, to_state, probability) -> None: ...
    def set_transition_pseudocount(self, from_state, to_state, count) -> None: ...
    def set_emission_score(self, seq_state, emission_state, probability) -> None: ...
    def set_emission_pseudocount(self, seq_state, emission_state, count) -> None: ...

class HiddenMarkovModel:
    state_alphabet: Incomplete
    emission_alphabet: Incomplete
    initial_prob: Incomplete
    transition_prob: Incomplete
    emission_prob: Incomplete
    def __init__(self, state_alphabet, emission_alphabet, initial_prob, transition_prob, emission_prob, transition_pseudo, emission_pseudo) -> None: ...
    def get_blank_transitions(self): ...
    def get_blank_emissions(self): ...
    def transitions_from(self, state_letter): ...
    def transitions_to(self, state_letter): ...
    def viterbi(self, sequence, state_alphabet): ...
