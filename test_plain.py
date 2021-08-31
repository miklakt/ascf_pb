# %%
import matplotlib.pyplot as plt
from ascf_pb.topology import kappa
from ascf_pb import factory
kwargs_ = dict(
    chi=0.3,
    N=1000,
    sigma = 0.02,
    R = 140
)

#kappa_cb = kappa.kappa_regular_dendron(2,4)
kappa_cb = kappa.linear
phi_profile = factory.phi(kappa_cb=kappa_cb, **kwargs_)
D = factory.D(kappa_cb=kappa_cb, **kwargs_)()
phi = [phi_profile(z) for z  in range(round(D+1.5))]
plt.plot(phi)
plt.show()