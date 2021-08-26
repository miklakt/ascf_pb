# %%
from ascf_pb.topology.pore import phi_D_unrestricted
import matplotlib.pyplot as plt
from ascf_pb.topology import kappa as Kappa
from ascf_pb.profiles import build_phi_profile_solver
from ascf_pb.topology import plain
chi=0.6
N=1000
sigma = 0.02
kappa_cb = Kappa.kappa_regular_dendron(2,4)
#kappa_cb = Kappa.kappa_plain
phi_profile = build_phi_profile_solver(
    kappa_cb, plain.D_universal, plain.phi_D_universal,
    chi = chi, N=N, sigma = sigma, R = 40)
phi = [phi_profile[0](z) for z  in range(round(phi_profile[1]+1.5))]
plt.plot(phi)
plt.show()