# %%
from ascf_pb.topology.kappa import kappa_plain
from ascf_pb.profiles import build_phi_profile_solver
from ascf_pb.topology.plain import D_universal, phi_D_universal
chi=0.6
N=1000
sigma = 0.02
kappa = kappa_plain(N)
phi_profile = build_phi_profile_solver(
    kappa_plain, D_universal, phi_D_universal,
    chi = chi, N=N, sigma = sigma, R = 50)
# %%
phi = [phi_profile(z) for z  in range(250)]
# %%
import matplotlib.pyplot as plt
plt.plot(phi)
# %%
from ascf_pb.topology.pore import D_unrestricted, phi_D_unrestricted
chi=0.6
N=1000
sigma = 0.02
Radius = 300
kappa = kappa_plain(N)
phi_profile = build_phi_profile_solver(
    kappa_plain, D_universal, phi_D_universal,
    chi = chi, N=N, sigma = sigma, R = 50)