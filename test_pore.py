#%%
import numpy as np
import matplotlib.pyplot as plt
from ascf_pb.solver import Phi
from ascf_pb.topology import pore
from ascf_pb.topology import kappa
from ascf_pb.profiles import build_phi_profile_solver
import matplotlib.pyplot as plt
def get_pore_phi_profile(N, sigma, chi, pore_Radius):
    phi_profile = build_phi_profile_solver(
        kappa.kappa_plain, pore.D_universal, pore.phi_D_universal,
        chi = chi, N=N, sigma = sigma, pore_Radius = pore_Radius, R = pore_Radius)
    return phi_profile
def draw_profile(phi_profile, pore_Radius):
    phi = [phi_profile[0](z) for z  in range(pore_Radius)]
    plt.margins(0,0)
    plt.xlim([0, pore_Radius])
    plt.ylim([0, phi_profile[0](0)])
    plt.plot(phi)
def draw_pore(phi_profile, pore_Radius):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.patches import Rectangle
    import mpl_toolkits.mplot3d.art3d as art3d
    import numpy as np
    fig = plt.figure()
    ax = Axes3D(fig)
    rad = np.arange(pore_Radius)
    azm = np.linspace(0, 2 * np.pi, 100)
    r, th = np.meshgrid(rad, azm)
    z = [[phi_profile[0](pore_Radius-i) for i in j] for j in r]

    polar_ax = plt.subplot(projection="polar")
    polar_ax.pcolormesh(th, r, z, rasterized = True, vmin = 0, vmax=0.3)

    polar_ax.set_xticks([])
    polar_ax.set_yticks([])

    channel_R = pore_Radius - phi_profile[1]
    if channel_R>0:
        polar_ax.plot(azm,[channel_R]*len(azm), color = 'red')

    return fig
#%%
chi = 0.45
sigma = 0.02
N=1000
pore_Radius = 150
phi_profile = get_pore_phi_profile(N,sigma, chi, pore_Radius)
draw_profile(phi_profile,pore_Radius)
# %%
chi = 0.7
phi_profile = get_pore_phi_profile(N,sigma, chi, pore_Radius)
fig = draw_pore(phi_profile, pore_Radius)
fig.savefig(f'pore_figs/N_{N}_sigma_{sigma}_chi_{chi}_Radius_{pore_Radius}.pdf')
# %%
