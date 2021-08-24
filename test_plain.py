# %%
import matplotlib.pyplot as plt
from ascf_pb.topology import kappa as Kappa
from ascf_pb.profiles import build_phi_profile_solver
from ascf_pb.topology import plain
chi=0.6
N=1000
sigma = 0.02
#kappa = Kappa.kappa_plain(N)
kappa = Kappa.kappa_regular_dendron(2,3)(N)
phi_profile = build_phi_profile_solver(
    Kappa.kappa_regular_dendron(1,3), plain.D_universal, plain.phi_D_universal,
    chi = chi, N=N, sigma = sigma, R = 150)
phi = [phi_profile[0](z) for z  in range(round(phi_profile[1]+1.5))]
plt.plot(phi)
plt.show()
# %%
norm = plain.normalization_unrestricted(chi, kappa, N*sigma)
# %%
plt.plot([norm(d) for d in range(100)])
# %%
int(phi_profile[1])*2
# %%
