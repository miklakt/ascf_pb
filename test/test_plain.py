# %%
import matplotlib.pyplot as plt
import numpy as np
import ascf_pb
import ascf_pb.topology
from functools import partial

brush = ascf_pb.BrushSolver()

chi=0.6
N=1000
sigma = 0.02
R = None

def get_by_kwargs(dataframe, **kwargs):
    return dataframe.loc[(dataframe[list(kwargs)] == pd.Series(kwargs)).all(axis=1)]

import pandas as pd
sfb_data = pd.read_pickle("empty_brush.pkl")
sfb_data = get_by_kwargs(sfb_data, chi_PS = chi, N=N)
#%%
#eta = ascf_pb.topology.regular_dendron_eta(1, 1)
D = brush.D(chi=chi, N=N, sigma = sigma, R=R)
z = np.arange(0, D)
phi_profile = brush.phi(chi=chi, N=N, sigma = sigma, R=R, z = z)
z = np.append(z, z[-1])
phi_profile = np.append(phi_profile, 0)

plt.plot(phi_profile, label = "SS-SCF")
plt.plot(sfb_data.phi.squeeze(), label = "SF-SCF")
plt.xlabel("z")
plt.ylabel("$\phi$")
plt.legend()
plt.show()
# %%
