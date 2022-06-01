#%%
import matplotlib.pyplot as plt
import ascf_pb
import ascf_pb.particle_geometry
import numpy as np
import pandas as pd

def get_by_kwargs(dataframe, **kwargs):
    return dataframe.loc[(dataframe[list(kwargs)] == pd.Series(kwargs)).all(axis=1)]

from scipy.optimize import least_squares
def fit_a1a2(sf_scf_pc, sfscf_free_energy, sfscf_phi, X0 = None, **kwargs):
    brush = ascf_pb.BrushInsertionFreeEnergy(external_phi=sfscf_phi)
    def resid(X):
        fe_estimated = brush.total_free_energy(**kwargs, pc = sf_scf_pc, expansion_coefs = tuple(X))
        r = fe_estimated - sfscf_free_energy
        return r
    if X0 is None:
        X0 = np.random.random(2)
    res = least_squares(resid, x0 =X0)
    return tuple(res.x)

chi=0.0
N=1000
sigma = 0.02
R = None
chi_PC =-1.0
#expansion_coefs = (0.19814812, -0.08488959)
ph = 4
pw = 4

# load sfbox data
sfscf_data = pd.read_pickle("empty_brush.pkl")
sfscf_phi = get_by_kwargs(sfscf_data, chi_PS = chi, N=N).phi.squeeze()
sfb_data = pd.read_pickle("post_proc.pkl")
plot_data = get_by_kwargs(sfb_data, chi_PS = chi, N=N, chi_PC =chi_PC, h=ph, w = pw)

brush = ascf_pb.BrushInsertionFreeEnergy()
brush_sfscf = ascf_pb.BrushInsertionFreeEnergy(external_phi=sfscf_phi)

fit_data = get_by_kwargs(sfb_data, chi_PS = chi, N=N, chi_PC =chi_PC, h=ph, w = pw)
fit_data=fit_data.loc[fit_data.phi>0.01]
sf_scf_pc= fit_data.ypos.values
sfscf_free_energy= fit_data.free_energy.values
expansion_coefs = fit_a1a2(
    fit_data.ypos.values, 
    fit_data.free_energy.values, 
    sfscf_phi, 
    pw = pw,
    ph = ph,
    chi=chi,
    chi_PC =chi_PC
    )


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
plt.suptitle("Precise and approximate free energy")

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
plt.suptitle("SS-SCF and SF-SCF free energy")


# %%
#plt.plot(z,PiV, label = "$\Pi$")
fig, ax = plt.subplots()
plt.plot(pc, tot_fe_sf, label = "$F(\phi_{SF-SCF})$")
plt.plot(pc, tot_fe, label = "$F(\phi_{SS-SCF})$", linestyle = ":")
#plt.plot(fapprx, label = "apprx")
plt.scatter(plot_data.ypos, plot_data.free_energy, label = "$F$ SF-SCF")
#plt.suptitle(f"{kwargs_} chi_PC: {chi_PC}")
plt.xlim(0, D*1.5)
plt.legend()
plt.suptitle("SS-SCF and SF-SCF(exact) free energy")
plt.text(  # position text relative to Axes
    0.05, 0.05, 'a1 = {:.2f}; a1 = {:.2f} '.format(*expansion_coefs),
    ha='left', va='bottom',
    transform=ax.transAxes
)
plt.show()

# %%
