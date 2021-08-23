<<<<<<< HEAD
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
    return d**2 - z**2 - (2/3)/kappa**2*(mu(phi_D,chi) + mu(phi,chi))

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
    def fsol(phi_ : float):
        return Z(
            phi = phi_, z = 0,
            d = d, kappa = kappa,
            chi = chi, phi_D = phi_D
        )
    return brentq(fsol, 0.99999, phi_D)

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

def phi_D_unrestricted(chi: float) -> float:
    """Calculates polymer volume fraction at the end of the brush`
    Args:
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: volume fraction
    """
    if chi <= 0.5:  # good solvent
        phi_D = 0
    else:  # poor solvent
        # to get phi_H we have to find where osmotic pressure vanishes
        def fsol(_phi):
            return Pi(_phi, chi)
        # find the root with brentq method
        min_phi = 0.00001  # exclude 0 from the roots
        max_phi = 0.99999  # exclude 1 from the roots
        try:
            phi_D = brentq(fsol, min_phi, max_phi)
        except Exception as e:
            raise e
    return phi_D

def normalization_plain_unrestricted(chi : float, kappa : float, theta : float):

    phi_D  = phi_D_unrestricted(chi)
    def integrand(z, d):
        return Phi(z = z, d=d,
        chi = chi, kappa=kappa, phi_D=phi_D)
    def integral(d):
        return integrate.quad(integrand, 0, d, args=(d,))[0] - theta
    return integral

def D_unrestricted(chi : float, kappa : float, theta : float, max_D_guess : float):
    normalization = normalization_plain_unrestricted(chi,kappa,theta)
    min_D = 0
    max_tries = 20
    tries = 0
    dD=max_D_guess/(max_tries+1)
    while tries < max_tries:
        try:

            _D = brentq(normalization, min_D, max_D_guess)
            break
        except ValueError as e:
            #print(e)
            #print(f'trying to decrease upper boundary max_H_guess')
            max_D_guess = max_D_guess-dD
            #print(f"max_H_guess : {max_H_guess}")
            tries = tries + 1
    else:
        raise ValueError()
    return _D
=======
from functools import lru_cache
from typing import Callable
import numpy as np
import scipy.integrate as integrate
from scipy import optimize

from .topology import kappa_plain

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

def Pi_phi_chi(phi: float, chi: float) -> float:
    """Calculate osmotic pressure for a given local polymer volume fraction
    Args:
        phi (float): local polymer volume fraction
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: osmotic pressure
    """
    Pi = -np.log(1-phi)-chi*phi**2-phi
    return Pi

def mu_phi_chi(phi : float, chi : float) -> float:
    """Chemical potential for a given volume fraction and solvent regime
    Args:
        phi (float): local polymer volume fraction
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: chemical potential
    """ 
    mu = -np.log(1-phi)-2*chi*phi
    return mu


def phi_H_free(chi: float) -> float:
    """Calculates polymer volume fraction at the end of the brush`
    Args:
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: volume fraction
    """
    if chi <= 0.5:  # good solvent
        phi_H = 0
    else:  # poor solvent
        # to get phi_H we have to find where osmotic pressure vanishes
        def fsol(_phi):
            return Pi_phi_chi(_phi, chi)
        # find the root with brentq method
        min_phi = 0.00001  # exclude 0 from the roots
        max_phi = 0.99999  # exclude 1 from the roots
        try:
            phi_H = optimize.brentq(fsol, min_phi, max_phi)
        except Exception as e:
            #print(f'optimize.brentq error in phi_H(chi = {chi})')
            raise e
    return phi_H

def phi_H_restricted(R: float) -> Callable[[float], float]:
    def phi_H(chi : float):
        _H = phi_H_free(chi)
        if _H <= R: 
            return _H
        else:
            return R



def _LAMBDA_2(H_2: float, chi: float, N: int, 
                K : Callable[[float], float], 
                phi_H : Callable[[float], float]
                ) -> float:
    """Calculates lambda_squared
    Args:
        H_2 (float): brush's height squared
        chi (float): Flory-Huggins parameter polymer-solvent
        N (int): polymer chain's length
    Returns:
        [type]: [description]
    """
    _phi_H = phi_H(chi)
    _mu = mu_phi_chi(_phi_H, chi)
    L_2 = H_2+ (2/3)*(1/K(N))**2 * _mu
    return L_2


def _Z_2(phi: float, chi: float, N: int, H_2: float, 
        K : Callable[[float], float], phi_H : Callable[[float], float]
        ) -> float:
    """Express z squared from strong-stretching SCF approximation
    Args:
        phi (float): local polymer volume fraction
        chi (float): Flory-Huggins parameter polymer-solvent
        N (int): polymer chain's length
        H_2 (float): brush's height squared
    Returns:
        (float): [description]
    """
    _mu = mu_phi_chi(phi, chi)
    L_2 = _LAMBDA_2(H_2, chi, N, K, phi_H)
    z_2 = L_2-(2/3)*(1/K(N))**2 * _mu
    return z_2


@lru_cache()
def _PHI_0(chi: float, N: int, H_2: float, 
            K : Callable[[float], float], phi_H : Callable[[float], float]
            ) -> float:
    """Calculates polymer volume fraction at the grafting surface
    Args:
        chi (float): Flory-Huggins parameter polymer-solvent
        N (int): polymer chain's length
        H_2 (float): brush's height squared
    Returns:
        float: distance from grafting surface squared
    """
    def fsol(_phi):
        return _Z_2(_phi, chi, N, H_2, K, phi_H)

    # the minimum value is at the brush's end
    min_phi = phi_H(chi)
    try:
        phi = optimize.brentq(fsol, min_phi, 0.99999)
    except Exception as e:
        #print(
        #    f'optimize.brentq error in _PHI_0(chi = {chi}, N = {N}, H_2={H_2})'
        #    )
        raise e
    return phi


def _Z_2_inv(z: float, chi: float, N: float, H: float, 
                K : Callable[[float], float], phi_H : Callable[[float], float]
            ) -> float:
    """Express phi from strong-stretching SCF approximation by inverting _Z_2
    Args:
        z (float): distance from grafting surface
        chi (float): Flory-Huggins parameter polymer-solvent
        N (float): polymer chain's length
        H_2 (float): brush's height squared
    Returns:
        float: local polymer volume fraction
    """
    z_2 = z**2
    H_2 = H**2
    if z_2 > H_2:  # outside the brush
        phi = 0
    elif z_2 == 0:  # at the grafting surface
        phi = _PHI_0(chi, N, H_2, K, phi_H)
    elif z_2 == H_2:  # at the brush's end
        phi = phi_H(chi)
    else:  # inside the brush
        # we can find correct phi by finding the root of fsol(_phi)
        def fsol(_phi):
            return _Z_2(_phi, chi, N, H_2, K, phi_H)-z_2

        max_phi = _PHI_0(chi, N, H_2, K, phi_H)
        min_phi = phi_H(chi)
        try:
            phi = optimize.brentq(fsol, min_phi, max_phi)
        except Exception as e:
            #print(
            #    f'optimize.brentq error in _Z_2_inv(z = {z}, chi = {chi}, N = {N}, H={H})'
            #)
            raise e
    return phi

@lru_cache()
def H(sigma: float, chi: float, N: int, 
        K : Callable[[float], float] = kappa_plain, 
        phi_H : Callable[[float], float] = phi_H_free,
        max_H_guess : float = None
        ) -> float:
    """Calculates brush's height (thickness) 
    by fullfiling normalization condition
    Args:
        sigma (float): grafting density (chains per square)
        chi (float): Flory-Huggins parameter polymer-solvent
        N (int): polymer chain's length
    Returns:
        float: brush's height
    """
    def fsol(z):  # normalization
        # integrate _Z_2_inv for a given sigma, chi, N from 0 to H
        return integrate.quad(_Z_2_inv, 0, z, args=(chi, N, z, K, phi_H))[0] - N*sigma
    min_H = 0.0
    #print("No max_H_guess provided")
    max_H_guess = N
    #print(f"max_H_guess : {max_H_guess}")
    max_tries = 20
    tries = 0
    dH=max_H_guess/(max_tries+1)
    while tries < max_tries:
        try:
            _H = optimize.brentq(fsol, min_H, max_H_guess)
            break
        except ValueError as e:
            #print(e)
            #print(f'trying to decrease upper boundary max_H_guess')
            max_H_guess = max_H_guess-dH
            #print(f"max_H_guess : {max_H_guess}")
            tries = tries + 1
    else:
        raise ValueError()
    return _H


def phi(sigma: float, chi: float, N: float, 
        K : Callable[[float], float] = kappa_plain,
        phi_H : Callable[[float], float] = phi_H_free
        ):
    """Calculates polymer volume fraction in a polymer brush at a given distance
    from the grafting surface
    Args:
        sigma (float): grafting density (chains per square)
        N (int): polymer chain's length
        chi (float): Flory-Huggins parameter polymer-solvent
        z (float): distance from the grafting surface
    Returns:
        function: volume fraction phi(z) function
    """

    _H = H(sigma, chi, N, K, phi_H)
    @np.vectorize
    def _phi(z : float) -> float:
        """volume fraction phi(z) function for given sigma, N, chi
        function is compatible with numpy.array input
        Args:
            z (float): distance from the grafting surface
        Returns:
            float: volume fraction
        """
        return _Z_2_inv(z, chi, N, _H, K, phi_H)
    return _phi
>>>>>>> c7e0c7d2be1708f9bb0529403ffa31666630f5a0
