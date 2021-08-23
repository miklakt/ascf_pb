#%%
from functools import lru_cache
from typing import Callable, Tuple
import numpy as np
import scipy
import scipy.integrate as integrate
from scipy import optimize

#from .topology import kappa_plain

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
    return d**2 - z**2 - (2/3)/kappa**2*(mu(phi_D,chi) + mu(phi,chi))

def phi_0(chi : float, kappa : float, d : float, phi_D : float):
    def fsol(phi_ : float):
        return Z(
            phi = phi_, z = 0,
            d = d, kappa = kappa,
            chi = chi, phi_D = phi_D
        )
    return optimize.brentq(fsol, 0.99999, phi_D)

def F(z : float, chi : float, kappa : float, d : float, phi_D : float):
    if z>d: return 0
    a = phi_0(chi, kappa, d, phi_D)
    if z==0: return a
    b = phi_D
    if z==phi_D: return b

    def fsol(phi_ : float):
        return Z(
            phi = phi_, z = z,
            d = d, kappa = kappa,
            chi = chi, phi_D = phi_D
        )
    return optimize.brentq(fsol, a, b)

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
            phi_D = optimize.brentq(fsol, min_phi, max_phi)
        except Exception as e:
            raise e
    return phi_D

def normalization_plain_unrestricted(chi : float, kappa : float, theta : float):
    phi_D  = phi_D_unrestricted(chi)
    def integrand(z, d):
        return F(z = z, d=d,
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

            _D = optimize.brentq(normalization, min_D, max_D_guess)
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


#%%
N =1000
chi =0
sigma = 0.02
theta = N*sigma
kappa = np.pi/(2*N)
#%%
D_unrestricted(chi,kappa,theta, N)




















# %%
F(100, 0.5, np.pi/(2*1000), 100,0)
# %%
normalization(0.5, np.pi/(2*1000), 0, 1000*0.02)(100.4)
# %%
F(101, 0, np.pi/(2*1000), 100, 0)

# %%
