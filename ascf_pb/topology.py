#%%
import numpy as np
import math
from typing import Tuple

def kappa_plain(N : float) -> float:
    k = 3/8 * np.pi**2/N**2
    return k

def kappa_regular_dendron_gnq(g : int, n : int, q :int) -> float:
    if g == 1:
        k = 1/n * math.atan(1/math.sqrt(q))
    if g == 2:
        k = 1/n * math.atan(1/math.sqrt(q*(q+2)))
    return k

def get_n(N : float, g : int, q : int) -> Tuple[float, float]:
    chains = 1+np.sum([q**i for i in range(g)]) #slow implementation
    n=N/chains
    return N

def kappa_regular_dendron(g : int, q : int):
    def kappa(N):
        n = get_n(N, g, q)
        return kappa_regular_dendron_gnq(g, n, q)
    return kappa