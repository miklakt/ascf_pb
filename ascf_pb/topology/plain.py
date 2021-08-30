from functools import lru_cache
from ascf_pb.topology.pore import phi_D_unrestricted
from ascf_pb.solver import Phi
from ascf_pb.topology import common
from scipy.optimize import brentq
from scipy import integrate
import numpy as np

required_keys = ['N', 'sigma', 'chi', 'z']

phi_D_unrestricted = common.phi_D_unrestricted

def theta(N : float, sigma : float):
    return N*sigma

def normalization_unrestricted(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        phi_D : float
        ):    
    def integrand(z, d):
        return Phi(z = z, d=d,
        chi = chi, kappa=kappa, phi_D=phi_D)
    def integral(d):
        return integrate.quad(integrand, 0, d, args=(d,))[0] - theta(N, sigma)
    return integral


def D_unrestricted(chi : float, kappa : float, N : float, sigma : float):
    phi_D = phi_D_unrestricted(chi)
    normalization = normalization_unrestricted(chi, kappa, N, sigma, phi_D)
    min_D = theta(N, sigma)
    if phi_D == 0: max_D = N
    else: max_D = theta(N, sigma)/phi_D
    _D = common.normalization_find_root(normalization, min_D, max_D)
    return _D


################################################################################
def normalization_restricted(
    chi : float, kappa : float, 
    N : float, sigma : float, 
    R : float
    ):
    def integrand(z, phi_R):
        return Phi(z = z, d=R,
        chi = chi, kappa=kappa, phi_D=phi_R)
    def integral(phi_R):
        return integrate.quad(integrand, 0, R, args = (phi_R,))[0] - theta(N, sigma)
    return integral


def phi_D_restricted(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float
        ):
    phi_R_min = phi_D_unrestricted(chi)
    phi_R_max = 0.99 #can be evaluated with a separate routine
    norm = normalization_restricted(chi, kappa, N, sigma, R)
    phi_R = brentq(norm, phi_R_min, phi_R_max)
    return phi_R


################################################################################
@lru_cache()
def phi_D_universal(chi : float, kappa : float, N : float, sigma : float, R : float = None):
    if R is None:
        phi_D = phi_D_unrestricted(chi)
        return phi_D
    #else
    try:
        print('Restriction notified.')
        D = D_unrestricted(chi, kappa, N, sigma)
        if R>=D or R<=0:
            #brush is not restricted
            print('Not restricted, D<=R.')
            phi_D = phi_D_unrestricted(chi)
        else:
            #brush is restricted
            print('Restricted, D>=R.')
            phi_D = phi_D_restricted(chi, kappa, N, sigma, R)
    except:
        phi_D = phi_D_restricted(chi, kappa, N, sigma, R)
    return phi_D

@lru_cache()
def D_universal(chi : float, kappa : float, N : float, sigma : float, R : float = None):
    D = D_unrestricted(chi, kappa, N, sigma)
    if R is None: return D
    if D<=R: 
        return D
    else:
        return R