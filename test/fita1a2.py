#%%
import matplotlib.pyplot as plt
import ascf_pb
import ascf_pb.particle_geometry
import numpy as np
import pandas as pd
plt.rcParams["figure.dpi"] = 300

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
ph = 4
pw = 4

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

fit_data = get_by_kwargs(sfscf_data, chi_PS = chi, N=N, chi_PC =chi_PC, h=ph, w = pw)
fit_data=fit_data.loc[fit_data.phi>0.01]
sf_scf_pc= fit_data.ypos.values
sfscf_free_energy= fit_data.free_energy.values

#%%
expansion_coefs = fit_a1a2(
    fit_data.ypos.values, 
    fit_data.free_energy.values, 
    sfscf_empty_phi, 
    pw = pw,
    ph = ph,
    chi=chi,
    chi_PC =chi_PC
    )

print (expansion_coefs)
#%%