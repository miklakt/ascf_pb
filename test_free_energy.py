#%%
import matplotlib.pyplot as plt
import ascf_pb
import ascf_pb.particle_geometry
import numpy as np

import logging
import sys
logging.basicConfig(level = logging.INFO, stream=sys.stdout)

kwargs_ = dict(
    chi=0.2,
    N=1000,
    sigma = 0.02,
    R = 100
)

ascf_pb.set_config(vectorize = True)

Pi_func = np.vectorize(ascf_pb.Pi(**kwargs_))
D = ascf_pb.D(**kwargs_)()
#%%
z  = np.arange(0, D+0.2, 0.1)
ph = 4
pw = 4
pz = np.arange(ph, D+ph, 0.5)
V = ascf_pb.particle_geometry.cylinder.volume(ph,pw)
pz = np.arange(ph, D+ph, 0.5)
osmotic_func = ascf_pb.osmotic_free_energy(ph=ph, pw=pw, **kwargs_)
#%%
Pi = Pi_func(z)
osmotic = osmotic_func
# %%
plt.plot(z,Pi, label = "$\Pi$")
plt.plot(pz,Pi_int, label = "osmotic")
plt.legend()
plt.show()

# %%
