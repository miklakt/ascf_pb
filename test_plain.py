# %%
import matplotlib.pyplot as plt
from ascf_pb.topology import kappa
from ascf_pb import factory
kwargs_ = dict(
    chi=0.2,
    N=1000,
    sigma = 0.02,
    R = 100
)

kappa_cb = kappa.kappa
eta = kappa.regular_dendron_eta(2,3)
phi_profile = factory.phi(kappa_cb=kappa_cb, eta=eta, **kwargs_)
D = factory.D(kappa_cb=kappa_cb, **kwargs_)()
phi = [phi_profile(z) for z  in range(round(D+1.5))]
plt.plot(phi)
plt.show()
# %%
