# %%
import matplotlib.pyplot as plt
import numpy as np
import ascf_pb
import ascf_pb.topology
import matplotlib.transforms as transforms
#plt.rcParams["figure.dpi"] = 300
#plt.rcParams['figure.figsize'] = [9.0, 6.0]

chi=0.5
N=300
sigma = 0.017
R = None

def get_by_kwargs(dataframe, **kwargs):
    return dataframe.loc[(dataframe[list(kwargs)] == pd.Series(kwargs)).all(axis=1)]

import pandas as pd
sfb_data = pd.read_pickle("empty_brush.pkl")
sfb_data = get_by_kwargs(sfb_data, chi_PS = chi, N=N)

brush = ascf_pb.BrushSolver()
brush_sfscf = ascf_pb.BrushInsertionFreeEnergy(external_phi=sfb_data.phi.squeeze())
#%%
%matplotlib ipympl
#eta = ascf_pb.topology.regular_dendron_eta(1, 1)
D = brush.D(chi=chi, N=N, sigma = sigma, R=R)
#percentile = 0.95
#D_sfscf = brush_sfscf.D(percentile =0.percentile)
z = np.arange(0, D)
phi_profile = brush.phi(chi=chi, N=N, sigma = sigma, R=R, z = z)
z = np.append(z, z[-1])
phi_profile = np.append(phi_profile, 0)
fig, ax = plt.subplots()

plt.plot(phi_profile, label = "SS-SCF")
#plt.plot(sfb_data.phi.squeeze(), label = "SF-SCF")
#plt.axvline(x = D_sfscf, linestyle = ":")
trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
#plt.text(s = f"percentile \n{percentile}", x = D_sfscf, y = 0.5, transform=trans)
plt.xlabel("z")
plt.ylabel("$\phi$")
plt.legend()
plt.show()
# %%