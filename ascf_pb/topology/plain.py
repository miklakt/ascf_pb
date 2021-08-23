from solver import Pi, Phi
from scipy.optimize import brentq
from scipy import integrate
import numpy as np

def phi_D_unrestricted(chi: float, **_) -> float:
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

def normalization_unrestricted(chi : float, kappa : float, theta : float):
    phi_D  = phi_D_unrestricted(chi)
    def integrand(z, d):
        return Phi(z = z, d=d,
        chi = chi, kappa=kappa, phi_D=phi_D)
    def integral(d):
        return integrate.quad(integrand, 0, d, args=(d,))[0] - theta
    return integral


def D_unrestricted(chi : float, kappa : float, N : float, sigma : float, **_):
    theta = N*sigma
    normalization = normalization_unrestricted(chi,kappa,theta)
    min_D = 0
    max_D = N
    max_tries = 20
    tries = 0
    dD=max_D/(max_tries+1)
    while tries < max_tries:
        try:

            _D = brentq(normalization, min_D, max_D)
            break
        except ValueError as e:
            #print(e)
            #print(f'trying to decrease upper boundary max_H_guess')
            max_D = max_D-dD
            #print(f"max_H_guess : {max_H_guess}")
            tries = tries + 1
    else:
        raise ValueError()
    return _D