"""
Provides routines to integrate particles over surface or volume
for cylinder
"""

import numpy as np

required_keys = ['ph', 'pw', "pc"]

def get_integration_interval(ph, pc):
    return pc-ph/2, pc+ph/2

def surface_integrand(ph, pw):
    def integrand(z):
        if z < 0 or z > ph:
            return 0
        elif z == 0 or z == ph:
            return np.pi*pw**2/4
        else:
            return np.pi*pw
    return integrand

def volume_integrand(ph, pw):
    def integrand(z):
        if z < 0 or z > ph:
            return 0
        else:
            return np.pi*pw**2/4
    return integrand

def volume(ph, pw):
    return np.pi*pw**2/4*ph

def surface(ph, pw):
    return np.pi*pw*ph+np.pi*pw**2/2
