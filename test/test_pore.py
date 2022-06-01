#%%
import matplotlib.pyplot as plt
import ascf_pb
import numpy as np
def draw_profile(z, phi, pore_Radius):
    z = np.append(z, z[-1])
    phi = np.append(phi, 0)
    plt.plot(z,phi)
    #plt.margins(0,0)
    plt.xlim([0, pore_Radius])
    #plt.ylim([0, max(phi)])

def draw_pore(phi_z, pore_Radius):
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
    z = [[phi_z(pore_Radius-i) for i in j] for j in r]

    polar_ax = plt.subplot(projection="polar")
    polar_ax.pcolormesh(th, r, z, rasterized = True, vmin = 0, vmax=0.3)

    polar_ax.set_xticks([])
    polar_ax.set_yticks([])

    channel_R = pore_Radius - D
    if channel_R>0:
        polar_ax.plot(azm,[channel_R]*len(azm), color = 'red')

    return fig
#%%
import ascf_pb
chi = 0.9
sigma = 0.02
N=2500
pore_Radius = 200
eta = 1
#eta = kappa.regular_dendron_eta(2,3)
brush = ascf_pb.BrushSolver("pore")
#%%
D = brush.D(N=N, 
    sigma=sigma, 
    chi=chi, 
    pore_Radius=pore_Radius, 
    R = pore_Radius)
z = np.linspace(0, D, 200)
phi= brush.phi(
        N=N, sigma=sigma, 
        chi=chi, 
        pore_Radius=pore_Radius, 
        R = pore_Radius, 
        eta = eta,
        z = z)
draw_profile(z, phi, pore_Radius)
#%%
phi_z = lambda z: brush.phi(
        N=N, sigma=sigma, 
        chi=chi, 
        pore_Radius=pore_Radius, 
        R = pore_Radius, 
        eta = eta,
        z = z)
fig = draw_pore(phi_z, pore_Radius)
#fig.savefig(f'pore_figs/N_{N}_sigma_{sigma}_chi_{chi}_Radius_{pore_Radius}.pdf')
# %%
