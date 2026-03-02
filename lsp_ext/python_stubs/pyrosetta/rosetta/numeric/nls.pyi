class lm_control_struct:
    epsilon: float
    ftol: float
    gtol: float
    maxcall: int
    printflags: int
    scale_diag: int
    stepbound: float
    xtol: float
    def __init__(self) -> None: ...

class lm_status_struct:
    fnorm: float
    info: int
    nfev: int
    def __init__(self) -> None: ...

def lm_enorm() -> float: ...
def lm_printout_std(n_par: int, par: float, m_dat: int, data: capsule, fvec: float, printflags: int, iflag: int, iter: int, nfev: int) -> None: ...
