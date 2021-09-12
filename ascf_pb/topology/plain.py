"""
The script provides routines to calculate plain topology restricted 
and unrestrestricted polymer brushes. 
The topology is taken into into accout in the normalization criteria routine
where an integrand were defined.

--------------------------------------------------------------------------------
                    UNRESTRICTED POLYMER BRUSH
There are routines to calculate brush thickness (D) and polymer density at the 
brush's end (phi_D) for an unrestricted brush. 
To get the brush thickness a root finding algorithm 
is used on the next expression (LaTeX):
   \int_{0}^{d}{\Phi(z, d, \chi, \kappa, \phi_D) dz} = \theta               (1)
   where brush thickness - d is the root of the expression.

For an unrestricted brush the polymer density at the brush's end (phi_D) is
defined by a vanishing osmotic pressure:
    \Pi(\chi, \phi_D) = 0                                                   (2)
    where polymer density phi_D is the root of the expression.
    phi_D = 0 for chi<0.5, while phi_D>0 for chi>0.5
--------------------------------------------------------------------------------
                    RESTRICTED POLYMER BRUSH
We impose restriction on a polymer brush by defining R, where R is the distance
from grafting surface to an impenetrable plain. Thus brush thickness is found
not from the normalization criteria (1), but is equal to R.

Slightly different normalization criteria is used to find polymer density at the 
brush's end (phi_D)
    \int_{0}^{R}{Phi(z, d=R, \chi, \kappa, phi_D) dz} = \theta              (3)
    where polymer density phi_D is the root of the expression.
--------------------------------------------------------------------------------
                    POSSIBLY RESTRICTED POLYMER BRUSH
Finally, we might not know if an imposed restriction is actually acts on a
polymer brush. This is easy to check by calculating the thickness of
the corresponding unrestricted brush (D_unrestricted). 
Obviously, if D_unrestricted > R the brush is restricted, so the routines for 
restricted brushes is employed. Otherwise, the routines for unrestricted 
polymer brushes is used.
--------------------------------------------------------------------------------

The script is imported in some other pieces of the code by string 'plain', 
to finnish the calculation the next keys are required:
    REQUIRED KEYS:
        N (float): the chain length of a polymer brush
        chi (float): Flory-Huggins parameter polymer-solvent
        z (float): distance from grafting surface
        sigma (float): grafting density (chains per unit area)

For specific cases of dendrons and restricted brushes the next keys 
has to be provided:
    OPTIONAL KEYS:
        kappa (float): topological parameter
        R (float): distance to an impenetrable surface 

@author: Mikhail Laktionov
miklakt@gmail.com
"""

from functools import lru_cache
from ascf_pb.topology.pore import phi_D_unrestricted
from ascf_pb.solver import Phi
from ascf_pb.topology import common
from scipy.optimize import brentq
from scipy import integrate
import numpy as np

# useful when importing from other pieces of code
required_keys = ['N', 'sigma', 'chi', 'z']

#when phi_D is defined by vanishing osmotic pressure
phi_D_unrestricted = common.phi_D_unrestricted

# amount of monomers in the system
def theta(N : float, sigma : float):
    return N*sigma

# this normalization is used to get brush thickness, 
# by finding d that satisfies normalization = theta
def normalization_unrestricted(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        phi_D : float
        ):
    # integrand for plain cartesian coordinates    
    def integrand(z, d):
        return Phi(z = z, d=d,
        chi = chi, kappa=kappa, phi_D=phi_D)
    # function for brentq to find d
    # which is the root of this integral
    def integral(d):
        return integrate.quad(integrand, 0, d, args=(d,))[0] - theta(N, sigma)
    # note that we return closure
    return integral

# get unrestricted brush thickness
def D_unrestricted(chi : float, kappa : float, N : float, sigma : float):
    # polymer density at the brush's end
    phi_D = phi_D_unrestricted(chi)
    # get normalization function for the brentq solver
    normalization = normalization_unrestricted(chi, kappa, N, sigma, phi_D)
    # lower boundary for D
    min_D = theta(N, sigma)
    # upper boundary for D
    if phi_D == 0: max_D = N
    else: max_D = theta(N, sigma)/phi_D
    # brentq with several tries
    _D = common.normalization_find_root(normalization, min_D, max_D)
    return _D


################################################################################
# this normalization is used to get polymer density at the brush's end for 
# a restricted brush, by finding phi_D that satisfies normalization = theta
def normalization_restricted(
    chi : float, kappa : float, 
    N : float, sigma : float, 
    R : float
    ):
    # integrand for plain cartesian coordinates
    def integrand(z, phi_R):
        return Phi(z = z, d=R,
        chi = chi, kappa=kappa, phi_D=phi_R)
    # function for brentq to find phi_D(phi_R)
    # which is the root of this integral
    def integral(phi_R):
        return integrate.quad(integrand, 0, R, args = (phi_R,))[0] - theta(N, sigma)
    return integral

# polymer density at the end of restricted brush
def phi_D_restricted(
        chi : float, kappa : float, 
        N : float, sigma : float, 
        R : float
        ):
    # lower boundary for phi_D, phi_D_restricted always bigger than phi_D_unrestricted
    phi_R_min = phi_D_unrestricted(chi)
    # upper boundary for phi_D
    phi_R_max = 0.99 #can be evaluated with a separate routine
    # get normalization function for the brentq solver
    norm = normalization_restricted(chi, kappa, N, sigma, R)
    # one try brentq,
    # there is several tries routine
    # but this case should not cause problems 
    phi_R = brentq(norm, phi_R_min, phi_R_max)
    return phi_R


################################################################################
# polymer density at the end of a brush, if we don't know if it is restricted
@lru_cache()
def phi_D_universal(chi : float, kappa : float, N : float, sigma : float, R : float = None):
    # if we have not provided restriction, it is unrestricted
    if R is None:
        phi_D = phi_D_unrestricted(chi)
        return phi_D
    # else, it might be restricted
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

# brush thickness, if we don't know if it is restricted
@lru_cache()
def D_universal(chi : float, kappa : float, N : float, sigma : float, R : float = None):
    D = D_unrestricted(chi, kappa, N, sigma)
    if R is None: return D
    if D<=R: 
        return D
    else:
        return R