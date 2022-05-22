"""
Routines to calculate insertion free energy of a probe particle
"""
import scipy.integrate as integrate
from typing import Callable
import numpy as np

required_keys = [
    "chi_PC",
    "expansion_coefs"
]

def osmotic_free_energy(z0, z1, volume_integrand, Pi_cb):
    if not isinstance(volume_integrand, Callable):
        def integrand(z_):
            return Pi_cb(z0)*volume_integrand
    # if A_cb provided as a callable
    else:
        def integrand(z_):
            return Pi_cb(z0)*volume_integrand(z_-z0)

    fe = integrate.quad(integrand, z0, z1)[0]
    return fe


def gamma_phi(phi, chi, chi_PC, expansion_coefs : list):
    chi_crit = 6*np.log(5/6)
    chi_ads = chi_PC - chi*(1-phi)
    psi = sum([a*phi**i for i, a in enumerate(expansion_coefs, start=1)])
    gamma = (chi_ads - chi_crit)*psi
    return gamma


def surface_free_energy(z0, z1, surface_integrand, gamma_cb, phi_cb):
    def integrand(z_):
        return surface_integrand(z_-z0)*gamma_cb(phi_cb(z_))
    fe = integrate.quad(integrand, z0, z1)[0]
    return fe