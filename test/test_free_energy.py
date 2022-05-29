#%%
import matplotlib.pyplot as plt
import ascf_pb
import ascf_pb.particle_geometry
import numpy as np

#import logging
#import sys
#logging.basicConfig(level = logging.INFO, stream=sys.stdout)

chi=0.6
N=1000
sigma = 0.02
R = None

phi_func = ascf_pb.phi(N=N, sigma=sigma, chi=chi, R=R)
Pi_func = ascf_pb.Pi(N=N, sigma=sigma, chi=chi, R=R)
D = ascf_pb.D(N=N, sigma=sigma, chi=chi, R=R)()
#%%
chi_PC =-0
expansion_coefs = [0.16, -0.08]
z  = np.arange(0, D+0.2, 0.1)
ph = 6
pw = 6
pz = np.arange(ph, D+ph, 0.5)
V = ascf_pb.particle_geometry.cylinder.volume(ph,pw)
S = ascf_pb.particle_geometry.cylinder.surface(ph,pw)
pz = np.arange(ph, D+ph, 0.5)
osmotic_func = ascf_pb.osmotic_free_energy(ph=ph, pw=pw, N=N, sigma=sigma, chi=chi, R=R)
surface_func = ascf_pb.surface_free_energy(
    N=N, sigma=sigma, chi=chi, R=R, 
    ph=ph, 
    pw=pw, 
    chi_PC = chi_PC, 
    expansion_coefs = expansion_coefs
    )
tot_func = ascf_pb.total_free_energy(
    N=N, sigma=sigma, chi=chi, R=R, 
    ph=ph, 
    pw=pw, 
    chi_PC = chi_PC, 
    expansion_coefs = expansion_coefs
)
#%%
PiV = Pi_func(z)*V
osmotic = osmotic_func(pc = pz)
surface = surface_func(pc = pz)
F = osmotic+surface

def get_by_kwargs(dataframe, **kwargs):
    return dataframe.loc[(dataframe[list(kwargs)] == pd.Series(kwargs)).all(axis=1)]

import pandas as pd
sfb_data = pd.read_pickle("post_proc.pkl")
plot_data = get_by_kwargs(sfb_data, chi_PS = chi, N=N, chi_PC =chi_PC, h=ph, w = pw)

from ascf_pb.free_energy import total_free_energy_apprx
fapprx = [total_free_energy_apprx(
    z=z_, surface = S, volume=V, phi_cb=phi_func, 
    Pi_cb=Pi_func, chi=chi, chi_PC = chi_PC, expansion_coefs = expansion_coefs
    ) for z_ in z]

# %%
plt.plot(z,PiV, label = "$\Pi$")
plt.plot(pz, F, label = "F/V")
plt.plot(pz, osmotic, label = "osm")
plt.plot(pz, surface, label = "sur")
plt.plot(fapprx, label = "apprx")
plt.scatter(plot_data.ypos, plot_data.free_energy)
#plt.suptitle(f"{kwargs_} chi_PC: {chi_PC}")
plt.legend()
plt.show()
# %%

# %%
