# %%
import matplotlib.pyplot as plt
import numpy as np
import ascf_pb
import ascf_pb.topology
from functools import partial

brush = ascf_pb.BrushSolver()

kwargs_ = dict(
    chi=0.2,
    N=1000,
    sigma = 0.02,
    R = 100
)
#%%
eta = ascf_pb.topology.regular_dendron_eta(2, 3)
D = brush.D(**kwargs_,eta=eta)
phi_profile = brush.phi(**kwargs_, z = np.linspace(0, D, 100), eta = eta)

plt.plot(phi_profile)
plt.show()
# %%
