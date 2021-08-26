from functools import lru_cache, partial
from typing import Callable, Tuple
import numpy as np
from scipy import integrate
from scipy.optimize import brentq
"""
The script is created to calculate volume fraction 
and osmotic pressure profiles of non-charged polymer brushes
using analytic self-consistent field (ASCF) method.
The parameters that describe the brush are the next ones:
chi -  Flory-Huggins parameter polymer-solvent
N - polymer chain's length
sigma - grafting density (chains per square)
z - distance from the grafting surface
The goal is functions that calculate volume fraction and osmotic pressure
with given parameters.
@author: Mikhail Laktionov
miklakt@gmail.com
"""

@lru_cache()
def Pi(phi: float, chi: float) -> float:
    """Calculate osmotic pressure for a given local polymer volume fraction
    Args:
        phi (float): local polymer volume fraction
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: osmotic pressure
    """
    Pi = -np.log(1-phi)-chi*phi**2-phi
    return Pi

@lru_cache()
def mu(phi : float, chi : float) -> float:
    """Chemical potential for a given volume fraction and solvent regime
    Args:
        phi (float): local polymer volume fraction
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: chemical potential
    """ 
    mu = -np.log(1-phi)-2*chi*phi
    return mu


def Z(phi : float,  d : float, z : float,  kappa : float, chi : float, phi_D : float):
    """An auxiliary function to calculate a polymer brush density profile.
    Profiles are founded by root finding of this function with all the arguments
    fixed, but phi

    Args:
        phi (float): polymer brush density
        d (float): polymer brush height
        z (float): distance from grafting surface
        kappa (float): topological parameter
        chi (float): Flory-Huggins parameter for the solvent-polymer interaction 
        phi_D (float): polymer density at the brush's end

    Returns:
        float: returns zero if arguments are consistent        
    """    
    return d**2 - z**2 + (2/3)/kappa**2*(mu(phi_D,chi) - mu(phi,chi))

@lru_cache()
def Phi_0(chi : float, kappa : float, d : float, phi_D : float):
    """Calculates polymer density at the grafting surface

    Args:
        chi (float): Flory-Huggins parameter for the solvent-polymer interaction
        kappa (float): topological parameter
        d (float): polymer brush height
        phi_D (float): polymer density at the brush's end
    Returns:
        float: polymer density at the grafting surface
    """
    eps=1e-09
    almost_one = 1-eps
    def fsol(phi_ : float):
        return Z(
            phi = phi_, z = 0,
            d = d, kappa = kappa,
            chi = chi, phi_D = phi_D
        )
    try:
        phi_0 = brentq(fsol, almost_one, phi_D)
    except ValueError as e:
        print("Warning brentq failed, phi must be too close to 1.0")
        phi_0 = almost_one
    return phi_0

@lru_cache()
def Phi(z : float, chi : float, kappa : float, d : float, phi_D : float):
    """Calculates polymer density at a given distance from the grafting surface
    Args:
        z (float): distance from grafting surface
        chi (float): Flory-Huggins parameter for the solvent-polymer interaction
        kappa (float): topological parameter
        d (float): polymer brush height
        phi_D (float): polymer density at the brush's end

    Returns:
        float: polymer density
    """    
    if z>d: return 0
    a = Phi_0(chi, kappa, d, phi_D)
    if z==0: return a
    b = phi_D
    if z==phi_D: return b

    def fsol(phi_ : float):
        return Z(
            phi = phi_, z = z,
            d = d, kappa = kappa,
            chi = chi, phi_D = phi_D
        )
    return brentq(fsol, a, b)


def min_mu(chi : float):
    if chi<=0.5:
        return 0
    else:
        return 1-np.log(1/(2*chi)) - 2*chi


def min_phi_D(chi : float):
    return max(0, 1-1/(2*chi))