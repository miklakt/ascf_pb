"""
Routines to calculate insertion free energy of a probe particle
"""
import scipy.integrate as integrate
from typing import Callable
import numpy as np
from functools import partial

required_keys = [
    "chi_PC",
    "expansion_coefs"
]


def osmotic_free_energy(z0 : float, z1 : float, volume_integrand : Callable, Pi_z : Callable) -> float:
    if not isinstance(volume_integrand, Callable):
        def integrand(z_ : float) -> float:
            return Pi_z(z0)*volume_integrand
    # if A_cb provided as a callable
    else:
        def integrand(z_ : float) -> float:
            return Pi_z(z0)*volume_integrand(z_-z0)
    fe = integrate.quad(integrand, z0, z1)[0]
    return fe


def gamma_phi(phi : float, chi : float, chi_PC : float, expansion_coefs : list) -> float:
    chi_crit = 6*np.log(5/6)
    chi_ads = chi_PC - chi*(1-phi)
    #psi = sum([a*phi**(i+1) for i, a in enumerate(expansion_coefs)])
    #gamma = (chi_ads - chi_crit)*psi
    a1, a2 = expansion_coefs
    psi = a1*phi + a2*phi**2
    gamma = (chi_ads-chi_crit)*psi
    return gamma


def surface_free_energy(
        z0 : float, z1 : float, 
        surface_integrand : Callable,
        A0 : float,
        A1 : float,
        phi_z : Callable, 
        chi : float, 
        chi_PC : float, 
        expansion_coefs : tuple
    ) -> float:
    gamma = partial(gamma_phi, 
        chi = chi, 
        chi_PC = chi_PC, 
        expansion_coefs = expansion_coefs
    )
    def integrand(z_):
        return surface_integrand(z_-z0)*gamma(phi_z(z_))
    fe = integrate.quad(integrand, z0, z1)[0] + A0*gamma(phi_z(z0)) + A1*gamma(phi_z(z1))
    return fe

def surface_free_energy_apprx(
    z: float,
    surface : float,
    phi_z : Callable,
    chi : float,
    chi_PC : float,
    expansion_coefs : tuple
)->float:
    phi = phi_z(z)
    gamma = gamma_phi(phi, chi, chi_PC, expansion_coefs)
    fe = surface*gamma
    return fe

def total_free_energy_apprx(
    z: float,
    surface : float,
    volume : float,
    phi_z : Callable,
    Pi_cb : Callable,
    chi : float,
    chi_PC : float,
    expansion_coefs : tuple
) -> float:
    phi = phi_z(z)
    Pi = Pi_cb(z)
    gamma = gamma_phi(phi, chi, chi_PC, expansion_coefs)
    fe = volume*Pi + surface*gamma
    return fe

def total_free_energy(
    z0 : float, 
    z1 : float, 
    surface_integrand : Callable, 
    volume_integrand : Callable, 
    phi_z : Callable, 
    Pi_cb : Callable,
    chi : float, 
    chi_PC : float,
    expansion_coefs : tuple
) -> float:
    osm = osmotic_free_energy(
        z0, z1, 
        volume_integrand, 
        Pi_cb
    )
    surf = surface_free_energy(
        z0, z1, 
        surface_integrand, 
        phi_z, 
        chi, 
        chi_PC, 
        expansion_coefs
    )

    return osm+surf