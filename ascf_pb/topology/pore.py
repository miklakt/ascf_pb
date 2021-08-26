from ascf_pb.solver import Pi, Phi
from ascf_pb.topology import utils
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


def normalization_unrestricted(chi : float, kappa : float, theta : float, pore_Radius : float):
    phi_D  = phi_D_unrestricted(chi)
    def integrand(z, d):
        return Phi(z = z, d=d,
        chi = chi, kappa=kappa, phi_D=phi_D)*abs(pore_Radius - z)
    def integral(d):
        return 2*np.pi*integrate.quad(integrand, 0, d, args=(d,))[0] - theta
    return integral


def D_unrestricted(
        chi : float, kappa : float,
        N : float, sigma : float, 
        pore_Radius : float, brentq_max_D = None,
        **_):
    theta = N*sigma*2*np.pi*pore_Radius
    normalization = normalization_unrestricted(chi,kappa,theta, pore_Radius)
    min_D = 0
    if brentq_max_D is None: max_D = min(pore_Radius,N)
    else: max_D = brentq_max_D
    _D = utils.normalization_find_root(normalization, min_D, max_D)
    return _D


################################################################################
def normalization_restricted(
        chi : float, kappa : float, 
        theta : float, R : float, 
        pore_Radius : float
        ):
    def integrand(z, phi_R):
        return Phi(z = z, d=R,
        chi = chi, kappa=kappa, phi_D=phi_R)*abs(pore_Radius - z)
    def integral(phi_R):
        return 2*np.pi*integrate.quad(integrand, 0, R, args = (phi_R,))[0] - theta
    return integral


def phi_D_restricted(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float, pore_Radius : float, 
        **_):
    phi_R_min = phi_D_unrestricted(chi)
    phi_R_max = 0.99 #can be evaluated with a separate routine
    theta = N*sigma*2*np.pi*pore_Radius
    norm = normalization_restricted(chi, kappa, theta, R, pore_Radius)
    phi_R = brentq(norm, phi_R_min, phi_R_max)
    return phi_R


################################################################################
def phi_D_universal(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float, pore_Radius : float, **_):
    try:
        D = D_unrestricted(chi, kappa, N, sigma, pore_Radius)
        print(R,D)
        if R>=D:
            #brush is not restricted
            print("Unrestricted")
            phi_D = phi_D_unrestricted(chi)
        else:
            #brush is restricted
            print("Restricted")
            phi_D = phi_D_restricted(chi, kappa, N, sigma, R, pore_Radius)
    except:
        #brush is restricted
        print("Restricted")
        phi_D = phi_D_restricted(chi, kappa, N, sigma, R, pore_Radius)
    return phi_D


def D_universal(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float, pore_Radius : float, **_):
    try:
        D_ = D_unrestricted(chi, kappa, N, sigma, pore_Radius)
        if D_>R: D_=R
    except:
        D_ = R
    return D_