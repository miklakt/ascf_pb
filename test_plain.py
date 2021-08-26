# %%
from ascf_pb.topology.pore import phi_D_unrestricted
import matplotlib.pyplot as plt
from ascf_pb.topology import kappa as Kappa
from ascf_pb.profiles import build_phi_profile_solver
from ascf_pb.topology import plain
chi=0.0
N=1000
sigma = 0.02
phi_profile = build_phi_profile_solver(
    Kappa.kappa_regular_dendron(2,3), plain.D_universal, plain.phi_D_universal,
    chi = chi, N=N, sigma = sigma, R=100)
phi = [phi_profile[0](z) for z  in range(round(phi_profile[1]+1.5))]
plt.plot(phi)
plt.show()
# %%
