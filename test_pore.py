#%%
import numpy as np
import matplotlib.pyplot as plt
from ascf_pb import profile_factory, topology
def draw_profile(phi_profile, pore_Radius):
    phi = [phi_profile(z) for z  in range(pore_Radius)]
    plt.margins(0,0)
    plt.xlim([0, pore_Radius])
    plt.ylim([0, phi_profile(0)])
    plt.plot(phi)
def draw_pore(phi_profile, pore_Radius, D):
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
    z = [[phi_profile(pore_Radius-i) for i in j] for j in r]

    polar_ax = plt.subplot(projection="polar")
    polar_ax.pcolormesh(th, r, z, rasterized = True, vmin = 0, vmax=0.3)

    polar_ax.set_xticks([])
    polar_ax.set_yticks([])

    channel_R = pore_Radius - D
    if channel_R>0:
        polar_ax.plot(azm,[channel_R]*len(azm), color = 'red')

    return fig
#%%
chi = 0.0
sigma = 0.02
N=1000
pore_Radius = 200
phi_profile = profile_factory.phi(topology='pore',N=N,sigma=sigma, chi=chi, pore_Radius=pore_Radius, R = pore_Radius)
D = profile_factory.D(topology='pore',N=N,sigma=sigma, chi=chi, pore_Radius=pore_Radius, R = pore_Radius)
draw_profile(phi_profile,pore_Radius)
fig = draw_pore(phi_profile, pore_Radius, D)
fig.savefig(f'pore_figs/N_{N}_sigma_{sigma}_chi_{chi}_Radius_{pore_Radius}.pdf')

# %%
