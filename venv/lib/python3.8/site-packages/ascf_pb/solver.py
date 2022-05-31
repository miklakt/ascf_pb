"""
Strong-Stretching Self Consistent Field formalism and core functions are 
defined in this file. The functions are mostly independent from 
the system geometry and a polymer brush topology.

Polymer density concentration or osmotic for a given parameters 
in the most cases can not be represented as an analytical function.

e.g. self-consistent molecular potential                        
- log(1 - φ) - 2⋅χ⋅φ = 3/2 κ^2 (Λ^2 - z^2)
φ(z) can not be expressed in analytical form here. Thus numerical root finding 
methods has to be exploit. 
The expression is implemented in the function Z(phi,  d, z, kappa, chi, phi_D),
solving Z=0 for a given argument gives us φ(z).


The main function script provides is Phi(z, chi, kappa, d, phi_D),
which solves aforesaid expression numerically for a given arguments.

Arguments used:
phi (float): local polymer volume fraction
chi (float): Flory-Huggins parameter polymer-solvent
d (float): polymer brush height
z (float): distance from grafting surface
kappa (float): topological parameter
phi_D (float): polymer density at the brush's end

@author: Mikhail Laktionov
miklakt@gmail.com
"""

from functools import lru_cache
import numpy as np
from scipy.optimize import brentq

class SolverError(Exception):
    """Custom error to help other components 
    to know where the raised error come from
    """    
    pass

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


def Z(
    phi : float,  d : float, z : float,  
    kappa : float, chi : float, phi_D : float
    ) -> float:
    """An auxiliary function to calculate a polymer brush density profile.
    Profiles are founded by root finding of this function with all the arguments
    fixed, but phi. See functions 'Phi_0' and 'Phi'.

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
def Phi_0(
    chi : float, kappa : float, d : float, 
    phi_D : float
    ) -> float:
    """Polymer density at the grafting surface. See function 'Z'

    Args:
        chi (float): Flory-Huggins parameter for the solvent-polymer interaction
        kappa (float): topological parameter
        d (float): polymer brush height
        phi_D (float): polymer density at the brush's end
    Returns:
        float: polymer density
    """
    #values close to 0 ant 1 to define open interval
    eps=1e-06
    almost_one = 1-eps
    
    #define function based on Z by set z=0 and fixing all the args, but phi
    def fsol(phi_ : float):
        return Z(
            phi = phi_, z = 0,
            d = d, kappa = kappa,
            chi = chi, phi_D = phi_D
        )
    try:
        #find f_0 that satisfies fsol(phi_0)=0
        phi_0 = brentq(fsol, almost_one, phi_D)
    except ValueError as e:
        print("Warning brentq failed, phi must be too close to 1.0")
        phi_0 = almost_one
    return phi_0

@lru_cache()
def Phi(
    z : float, chi : float, kappa : float, 
    d : float, phi_D : float) -> float:
    """Polymer density at a given distance from the grafting surface.
    See function 'Z'
    Args:
        z (float): distance from grafting surface
        chi (float): Flory-Huggins parameter for the solvent-polymer interaction
        kappa (float): topological parameter
        d (float): polymer brush height
        phi_D (float): polymer density at the brush's end

    Returns:
        float: polymer density
    """
    #if outside brush, phi=0    
    if z>d: return 0
    #phi(z=0) - density at the grafting surface 
    a = Phi_0(chi, kappa, d, phi_D)
    if z==0: return a
    #phi(z=D) - density at the brush's end
    b = phi_D
    if z==phi_D: return b

    #define function based on Z by fixing all the args, but phi
    def fsol(phi_ : float):
        return Z(
            phi = phi_, z = z,
            d = d, kappa = kappa,
            chi = chi, phi_D = phi_D
        )
    try:
        #find f_0 that satisfies fsol(phi_0)=0
        #we also know that phi(z=0)>phi(z)>phi(z=D)
        return brentq(fsol, a, b)
    except ValueError as e:
        #args must hav been incompatible, e.g. phi_D is to high
        print("solver.Phi() brentq failed")
        raise SolverError("solver.Phi() brentq failed")


def min_mu(chi : float) -> float:
    """Minimal molecular(segment) potential possible, 
    minimal value of -log(1-φ)-2⋅χ⋅φ 

    Args:
        chi (float): Flory-Huggins parameter for the solvent-polymer interaction

    Returns:
        float: molecular(segment) potential
    """    
    if chi<=0.5:
        return 0
    else:
        return 1-np.log(1/(2*chi)) - 2*chi


def min_phi_D(chi : float) -> float:
    """Polymer density corresponds to minimal molecular(segment) potential, 
    argmin(-log(1-φ)-2⋅χ⋅φ).

    Args:
        chi (float): Flory-Huggins parameter for the solvent-polymer interaction

    Returns:
        float: polymer brush density
    """
    return max(0, 1-1/(2*chi))