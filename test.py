# %%
from ascf_pb.topology.plain import phi_D_restricted, phi_D_unrestricted, D_unrestricted
from ascf_pb.topology.kappa import kappa_plain
from ascf_pb.solver import Phi_0
chi=0
N=1000
sigma = 0.02
kappa = kappa_plain(N)
#%%
phi_D = phi_D_unrestricted(chi)
D = D_unrestricted(chi, kappa, N, sigma)
phi_0 = Phi_0(chi, kappa, D, phi_D)
print(phi_0, phi_D)
#%%
R=D*0.5
phi_D_r = phi_D_restricted(chi, kappa, N, sigma, R=R)
phi_0_r = Phi_0(chi, kappa, R, phi_D_r)
print(phi_0_r, phi_D_r)
# %%
from ascf_pb.solver import Phi
from functools import partial
import numpy as np
phi = partial(Phi,d=D, chi = chi, kappa=kappa, phi_D=phi_D)
z = np.arange(round(D+1.5))
phi_ = [phi(z_) for z_ in z]
import matplotlib.pyplot as plt
plt.plot(z, phi_)

phi = partial(Phi,d=R, chi = chi, kappa=kappa, phi_D=phi_D_r)
z = np.arange(round(R+1.5))
phi_ = [phi(z_) for z_ in z]
plt.plot(z, phi_)
plt.show()
