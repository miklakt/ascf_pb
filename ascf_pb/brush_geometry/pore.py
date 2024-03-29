"""
The script provides routines to calculate pore topology restricted 
and unrestrestricted polymer brushes. 

Docs for this module are not yet written, please take a look on 
plain.py script in the same directory, to get the idea of what is going on

@author: Mikhail Laktionov
miklakt@gmail.com
"""

from functools import lru_cache

from ascf_pb.solver import Phi
from . import common
from scipy.optimize import brentq
from scipy import integrate
import numpy as np

phi_D_unrestricted = common.phi_D_unrestricted

def theta(N : float, sigma : float, pore_Radius : float):
    return N*sigma*2*np.pi*pore_Radius

def normalization_unrestricted(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        phi_D : float, pore_Radius : float):
    def integrand(z, d):
        return Phi(z = z, d=d,
        chi = chi, kappa=kappa, phi_D=phi_D)*abs(pore_Radius - z)
    def integral(d):
        return 2*np.pi*integrate.quad(integrand, 0, d, args=(d,))[0] - theta(N, sigma, pore_Radius)
    return integral

def D_boundary(pore_Radius : float, phi_const : float, N : float, sigma : float):
    D = pore_Radius - np.sqrt(pore_Radius**2 - theta(N, sigma, pore_Radius)/(np.pi*phi_const))
    return D

def D_unrestricted(
        chi : float, kappa : float,
        N : float, sigma : float, 
        pore_Radius : float):
    phi_D  = phi_D_unrestricted(chi)
    normalization = normalization_unrestricted(
                        chi, kappa, N, sigma, phi_D, pore_Radius)
    min_D = D_boundary(pore_Radius, 1, N, sigma)
    #min_D = 0
    if phi_D == 0: max_D=min(pore_Radius,N)
    else: max_D=min(pore_Radius,D_boundary(pore_Radius, phi_D, N, sigma))#min(pore_Radius,N)#
    _D = common.normalization_find_root(normalization, min_D, max_D)
    return _D


################################################################################
def normalization_restricted(
        chi : float, kappa : float, 
        N : float, sigma : float, R : float, 
        pore_Radius : float
        ):
    def integrand(z, phi_R):
        return Phi(z = z, d=R,
        chi = chi, kappa=kappa, phi_D=phi_R)*abs(pore_Radius - z)
    def integral(phi_R):
        return 2*np.pi*integrate.quad(integrand, 0, R, args = (phi_R,))[0] - theta(N, sigma, pore_Radius)
    return integral


def phi_D_restricted(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float, pore_Radius : float):
    phi_R_min = phi_D_unrestricted(chi)
    phi_R_max = 0.99 #can be evaluated with a separate routine
    norm = normalization_restricted(chi, kappa, N, sigma, R, pore_Radius)
    phi_R = brentq(norm, phi_R_min, phi_R_max)
    return phi_R


################################################################################
#@lru_cache()
def phi_D_universal(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float, pore_Radius : float):
    try:
        D = D_unrestricted(chi, kappa, N, sigma, pore_Radius)
        if R>=D:
            #brush is not restricted
            logging.debug("Unrestricted")
            phi_D = phi_D_unrestricted(chi)
        else:
            #brush is restricted
            logging.debug("Restricted")
            phi_D = phi_D_restricted(chi, kappa, N, sigma, R, pore_Radius)
    except:
        #brush is restricted
        logging.debug("Restricted")
        phi_D = phi_D_restricted(chi, kappa, N, sigma, R, pore_Radius)
    return phi_D

#@lru_cache()
def D_universal(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float, pore_Radius : float):
    try:
        D_ = D_unrestricted(chi, kappa, N, sigma, pore_Radius)
        if D_>R: D_=R
    except:
        D_ = R
    return D_


################################################################################
def normalization_pore_opening(
    chi : float, kappa : float, 
    N : float, sigma : float, phi_D : float
    ):
    def integrand(z, d):
        return Phi(z = z, d=d,
            chi = chi, kappa=kappa, phi_D=phi_D)*abs(d - z)
    def integral(d):
        return 2*np.pi*integrate.quad(integrand, 0, d, args=(d,))[0] - theta(N, sigma, d)
    return integral

#@lru_cache()
def opening_pore_Radius(
    chi : float, kappa : float,
    N : float, sigma : float,
    ):
    phi_D = phi_D_unrestricted(chi)
    min_pore_R = 2*N*sigma
    if phi_D == 0: max_pore_R=N
    else: max_pore_R=2*N*sigma/phi_D
    normalization = normalization_pore_opening(chi, kappa, N, sigma, phi_D)
    pore_R = common.normalization_find_root(normalization, min_pore_R, max_pore_R)
    return pore_R


def normalization_chi_opening(
    kappa : float, 
    N : float, sigma : float,
    R : float
    ):
    def integrand(z, chi):
        phi_D = phi_D_unrestricted(chi)
        return Phi(z = z, d=R,
            chi = chi, kappa=kappa, phi_D=phi_D)*abs(R - z)
    def integral(chi):
        return 2*np.pi*integrate.quad(integrand, 0, R, args=(chi,))[0] - theta(N, sigma, R)
    return integral

#@lru_cache()
def chi_opening(
    kappa : float,
    N : float, sigma : float, R : float
    ):
    chi_min=0
    chi_max=1
    normalization = normalization_chi_opening(kappa, N, sigma, R)
    try:
        chi_open = brentq(normalization, chi_min, chi_max)
    except:
        chi_open = None
    return chi_open
