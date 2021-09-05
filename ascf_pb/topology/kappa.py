#%%
import numpy as np
import math

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



###probably be removed later
def _dendron_gnq(g : int, n : int, q :int):
    if g == 1:
        k = 1/n * math.atan(1/math.sqrt(q))
    if g == 2:
        k = 1/n * math.atan(1/math.sqrt(q*(q+2)))
    return k

def _get_n(N : float, g : int, q : int):
    n_branches = 1+np.sum([q**i for i in range(1,g+1)]) #slow implementation
    n=N/n_branches
    return n

def _regular_dendron(g : int, q : int):
    def kappa(N):
        n = _get_n(N, g, q)
        return _dendron_gnq(g, n, q)
    return kappa