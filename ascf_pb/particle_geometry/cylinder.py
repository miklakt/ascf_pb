#%%
"""
Provides routines to integrate particles over surface or volume
for cylinder
"""

import numpy as np

required_keys = ['ph', 'pw', "pc"]

def integration_interval(ph : float, pc : float) -> tuple:
    return pc-ph/2, pc+ph/2

def surface_integrand(ph : float, pw : float):
    def integrand(z):
        if z < 0 or z > ph:
            return 0
        else:
            return np.pi*pw
    A0 = A1 = np.pi*pw**2/4
    return integrand, A0, A1

def volume_integrand(ph : float, pw : float):
    def integrand(z : float) -> float:
        if z < 0 or z > ph:
            return 0
        else:
            return np.pi*pw**2/4
    return integrand

def volume(ph : float, pw : float) -> float:
    return np.pi*pw**2/4*ph

def surface(ph : float, pw : float) -> float:
    return np.pi*pw*ph+np.pi*pw**2/2


#%%
if __name__ == "__main__":
    from scipy import integrate
    ph=4
    pw =4

    S = surface(ph, pw)
    integrand, A0, A1 = surface_integrand(ph, pw)
    S_int = integrate.quad(integrand, 0, ph)[0] + A0 + A1
# %%
