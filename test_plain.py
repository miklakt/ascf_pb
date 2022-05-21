# %%
import matplotlib.pyplot as plt
import numpy as np
from ascf_pb.topology import kappa
import ascf_pb

ascf_pb.set_config(vectorize = True)

kwargs_ = dict(
    chi=0.2,
    N=1000,
    sigma = 0.02,
    R = 100
)

kappa_cb = kappa.kappa
eta = kappa.regular_dendron_eta(2,3)
phi_profile = ascf_pb.phi(kappa_cb=kappa_cb, eta=eta, **kwargs_)
D = ascf_pb.D(kappa_cb=kappa_cb, **kwargs_)()
z = np.arange(0, D)
phi = phi_profile(z)
plt.plot(phi)
plt.show()
# %%
