# %%
import matplotlib.pyplot as plt
from ascf_pb.topology import kappa
from ascf_pb import profile_factory
kwargs_ = dict(
    chi=0.3,
    N=1000,
    sigma = 0.02,
    R = 140
)

#kappa_cb = kappa.kappa_regular_dendron(2,4)
kappa_cb = kappa.kappa_plain
phi_profile = profile_factory.phi(kappa_cb=kappa_cb, **kwargs_)
D = profile_factory.D(kappa_cb=kappa_cb, **kwargs_)()
phi = [phi_profile(z) for z  in range(round(D+1.5))]
plt.plot(phi)
plt.show()
