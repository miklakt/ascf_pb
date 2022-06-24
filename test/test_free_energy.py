#%%
import matplotlib.pyplot as plt
import ascf_pb
import ascf_pb.particle_geometry
import numpy as np
import pandas as pd
plt.rcParams["figure.dpi"] = 300
plt.rcParams['figure.figsize'] = [12.0, 8.0]

def get_by_kwargs(dataframe, **kwargs):
    return dataframe.loc[(dataframe[list(kwargs)] == pd.Series(kwargs)).all(axis=1)]

chi=0.6
N=1000
sigma = 0.02
R = None
chi_PC =-1.0
expansion_coefs = (0.198, -0.085)
ph = 10
pw = 10

# load sfbox data
# empty brush
sfscf_empty_data = pd.read_pickle("empty_brush.pkl")
sfscf_empty_phi = get_by_kwargs(sfscf_empty_data, chi_PS = chi, N=N).phi.squeeze()
# brush with an inserted particle
sfscf_data = pd.read_pickle("post_proc.pkl")
sfscf_data = get_by_kwargs(sfscf_data, chi_PS = chi, N=N, chi_PC =chi_PC, h=ph, w = pw)

# pure analytical SS-SCF brush
brush = ascf_pb.BrushInsertionFreeEnergy()
# analytical formulae applied to volume fraction profile calculated via SF_SCF
brush_sfscf = ascf_pb.BrushInsertionFreeEnergy(external_phi=sfscf_empty_phi)


D = brush.D(N=N, sigma=sigma, chi=chi, R=R)
z  = np.arange(0, D*1.5, 0.25)
pc = np.arange(ph, D*1.5, 0.25)

V = brush.Particle.volume(ph, pw)
S = brush.Particle.surface(ph, pw)

#%%
################################################################################
###Precise and approximate free energy comparison###############################
################################################################################
osmotic_fe = brush.osmotic_free_energy(
    ph=ph, pw=pw, N=N, sigma=sigma, chi=chi, R=R, pc = pc
    )
surface_fe = brush.surface_free_energy(
    N=N, sigma=sigma, chi=chi, R=R, 
    ph=ph, 
    pw=pw, 
    chi_PC = chi_PC, 
    expansion_coefs = expansion_coefs,
    pc = pc
    )
tot_fe = osmotic_fe+surface_fe
PiV = brush.Pi(N=N, sigma=sigma, chi=chi, R=R, z =z) * V

phi_z = lambda z_: brush.phi(N=N, sigma=sigma, chi=chi, R=R, z=z_)
surface_fe_apprx = [ascf_pb.free_energy.surface_free_energy_apprx(z_, S, phi_z, chi, chi_PC, expansion_coefs) for z_ in z]

plt.plot(z,PiV, label = "$\Pi V$", linestyle = ":")
plt.plot(pc, osmotic_fe, label = "osmotic")
plt.plot(z, surface_fe_apprx, label = "$\gamma S$", linestyle = ":")
plt.plot(pc, surface_fe, label = "surface")
plt.legend()
plt.suptitle(
    "Precise and approximate free energy\n" +\
    f"chi: {chi} chi_PC: {chi_PC} ph: {ph} pw: {pw}"
    )

#%%
################################################################################
###SS-SCF and SF-SCF free energy comparison#####################################
################################################################################
osmotic_fe_sf = brush_sfscf.osmotic_free_energy(
    chi = chi,
    ph=ph, pw=pw, pc = pc
    )
surface_fe_sf= brush_sfscf.surface_free_energy(
    chi = chi,
    chi_PC = chi_PC,
    ph=ph, 
    pw=pw,
    expansion_coefs = expansion_coefs,
    pc = pc
    )

tot_fe_sf = osmotic_fe_sf+surface_fe_sf

#plt.plot(z,PiV, label = "$\Pi V$")
plt.plot(pc, osmotic_fe, label = "osmotic SS-SCF", linestyle = ":")
plt.plot(pc, surface_fe, label = "surface SS-SCF", linestyle = ":")
plt.plot(pc, tot_fe, label = "total SS-SCF", linestyle = ":")

plt.plot(pc, osmotic_fe_sf, label = "osmotic SF-SCF")
plt.plot(pc, surface_fe_sf, label = "surface SF-SCF")
plt.plot(pc, tot_fe_sf, label = "total SF-SCF")
plt.axhline(y=0)
plt.legend()
plt.suptitle(
    "SS-SCF and SF-SCF free energy\n" +\
    f"chi: {chi} chi_PC: {chi_PC} ph: {ph} pw: {pw}"
    )


# %%
#plt.plot(z,PiV, label = "$\Pi$")
fig, ax = plt.subplots()
plt.plot(pc, tot_fe_sf, label = "$F(\phi_{SF-SCF})$")
plt.plot(pc, tot_fe, label = "$F(\phi_{SS-SCF})$", linestyle = ":")
#plt.plot(fapprx, label = "apprx")
plt.scatter(sfscf_data.ypos, sfscf_data.free_energy, label = "$F$ SF-SCF")
#plt.suptitle(f"{kwargs_} chi_PC: {chi_PC}")
plt.xlim(0, D*1.5)
plt.legend()
plt.suptitle(
    "SS-SCF and SF-SCF(exact) free energy\n" +\
    f"chi: {chi} chi_PC: {chi_PC} ph: {ph} pw: {pw}"
    )
plt.text(  # position text relative to Axes
    0.05, 0.05, 'a1 = {:.3f}; a1 = {:.3f} '.format(*expansion_coefs),
    ha='left', va='bottom',
    transform=ax.transAxes
)
plt.show()

# %%
