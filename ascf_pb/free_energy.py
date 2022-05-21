"""
Routines to calculate insertion free energy of a probe particle
"""

import scipy.integrate as integrate
from typing import Callable
from functools import reduce
import numpy as np

def osmotic_free_energy(z0, z1, A_cb, Pi_cb):
    if not isinstance(A_cb, Callable):
        def integrand(z_):
            return Pi_cb(z0)*A_cb
    #if A_cb provided as callable
    else:
        def integrand(z_):
            return Pi_cb(z0)*A_cb(z_-z0)
    
    fe = integrate.quad(integrand, z0, z1)[0]
    return fe

def gamma_phi(phi, chi_PS, chi_PC, a_virial):
    chi_crit = 6*np.log(5/6)
    chi_ads = chi_PC - chi_PS*(1-phi)
    psi = sum([a*phi**i for i, a in enumerate(a_virial, start = 1)])
    gamma = (chi_ads - chi_crit)*psi
    return gamma