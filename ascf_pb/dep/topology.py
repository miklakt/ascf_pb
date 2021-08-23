import numpy as np
import math
from typing import Tuple

def kappa_plain(N : float) -> float:
    k = np.pi/(2*N)
    return k

def kappa_regular_dendron_gnq(g : int, n : int, q :int) -> float:
    if g == 1:
        k = 1/n * math.atan(1/math.sqrt(q))
    if g == 2:
        k = 1/n * math.atan(1/math.sqrt(q*(q+2)))
    return k

def get_n(N : float, g : int, q : int) -> Tuple[float, float]:
    n_branches = 1+np.sum([q**i for i in range(1,g+1)]) #slow implementation
    n=N/n_branches
    return n

def kappa_regular_dendron(g : int, q : int):
    def kappa(N):
        n = get_n(N, g, q)
        return kappa_regular_dendron_gnq(g, n, q)
    return kappa