import numpy as np
import math

_arg_description = dict(
    eta="topological parameter",
    g = "dendron generations",
    q = "dendron functionality"
)


def kappa(N : float, eta : float = 1.0) -> float:
    kappa = np.pi*eta/(2*N)
    return kappa

def regular_dendron_eta(g : int, q : int) -> float:
    qq = 1+np.sum([q**i for i in range(1,g+1)])
    if g == 1:
        eta = math.atan(1/math.sqrt(q))*2/np.pi*qq
    if g == 2:
        eta =  math.atan(1/math.sqrt(q*(q+2)))*2/np.pi*qq
    return eta