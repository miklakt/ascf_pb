#%%
import matplotlib.pyplot as plt
import ascf_pb
from ascf_pb.topology import kappa
from ascf_pb import factory
import numpy as np
kwargs_ = dict(
    chi=0.2,
    N=1000,
    sigma = 0.02,
    R = 100
)

phi_profile = np.vectorize(factory.phi(**kwargs_))
Pi_profile = np.vectorize(factory.Pi(**kwargs_))
D = factory.D(**kwargs_)()
z  = np.arange(0, D+0.2, 0.1)
phi = phi_profile(z)
Pi = Pi_profile(z)

h = 4
A=16
V = A*h
pz = np.arange(h, D+h, 0.5)
Pi_int = np.array(
        [
        ascf_pb.osmotic_free_energy(z_-h/2, z_+h/2, A, Pi_profile) for z_ in pz
        ]
    )/V
# %%
plt.plot(z,Pi, label = "$\Pi$")
plt.plot(pz,Pi_int, label = "osmotic")
plt.legend()
plt.show()

# %%
