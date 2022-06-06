#%%
import ascf_pb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
#plt.rcParams["figure.dpi"] = 300
#plt.rcParams['figure.figsize'] = [12.0, 8.0]
LAST_USED_COLOR = lambda: plt.gca().lines[-1].get_color()

Gamma = 4.9 #pmol/cm2
D_experimental = 25 #nm
a = 0.76 #nm
N = 300
NA = 6.02214 #1e-23 mol-1
#%%
sigma = Gamma*a**2*NA*1e-3
# %%
brush = ascf_pb.BrushSolver()
from scipy import optimize
def find_chi_on_D(D, D_error, **kwargs):
    f = [
        lambda x: brush.D(chi = x, **kwargs)/a-D,
        lambda x: brush.D(chi = x, **kwargs)/a-D-D_error,
        lambda x: brush.D(chi = x, **kwargs)/a-D+D_error,
    ]
    chi = [optimize.brentq(f_, a=0.0, b=1.0) for f_ in f]
    return chi

chi = np.linspace(0, 1, 100)
D_on_chi = brush.D(N=N, sigma = sigma, chi = chi) / a
D_min = N*sigma/a
D_err = 2 #nm
d_chi = np.gradient(chi, D_on_chi)
chi_err = d_chi*D_err
chi_expected, chi_down, chi_up = find_chi_on_D(D_experimental, D_err, N=N, sigma=sigma)
phi_av = N*sigma/(D_experimental/a)
phi_err = N*sigma/(D_experimental/a)**2*D_err

rho = 1400 #mg/ml
#conc = phi_av*rho
#conc_err = phi_err*rho

M=64.1 #kg/mol
conc = Gamma*M/D_experimental*10 #mg/ml
# %%
%matplotlib ipympl
fig, ax = plt.subplots()
ax.plot(D_on_chi, chi)
ax.plot(D_on_chi, chi-chi_err, linestyle = ":", color = LAST_USED_COLOR())
ax.plot(D_on_chi, chi+chi_err, linestyle = ":", color = LAST_USED_COLOR())
ax.set_ylabel("$\chi$")
ax.set_xlabel("D, [nm]")
ax.set_xlim(0,max(D_on_chi))
ax.axvline(x = D_min, linestyle = ":", color = "black")
ax.axhline(y = chi_expected, color = "red", linewidth = 0.2)
ax.scatter(x = D_experimental, y=chi_expected, marker = "o", color = "red")
ax.vlines(x = D_experimental, color = "red", ymin = chi_down, ymax = chi_up)
ax.text(x = D_experimental, color = "red", y = chi_down, s = '{:.3f}'.format(chi_down), va = "top")
ax.text(x = D_experimental, color = "red", y = chi_up, s = '{:.3f}'.format(chi_up), va = "bottom")
trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
ax.text(s="$D_{min} \equiv N\sigma / a$", x = D_min, y=0, va = "bottom", ha = "left", transform = trans)
annotation_text = f"D = {D_experimental} $\pm$ {D_err} nm" +\
    f"\n$\Gamma = {Gamma}$ pmol" +\
    f"\n$a = {a}$ nm" +\
    f"\n$N = {N}$" +\
    f"\n$\sigma = {sigma :.4f}$" +\
    "\n$\overline{\phi} =$"+f"{phi_av :.2f}" + f"$\pm {phi_err :.2f}$" +\
    "\n$C =$"+f"{conc :.3f}"+ "mg/ml"
    #f"$C/\rho = {conc/rho :.3f}$"

ax.text(s=annotation_text, x = 0.98, y=0.98, va = "top", ha = "right", transform = ax.transAxes)
trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
ax.text(x = 0.98, color = "red", y = chi_expected, s = "$\chi =$"+'{:.3f}'.format(chi_expected), va = "bottom", ha = 'right', transform = trans)
# %%
fig, ax = plt.subplots()
ax.plot(D_on_chi, d_chi)
ax.set_ylabel(r"$\partial \chi / \partial D$")
ax.set_xlabel("D, [nm]")
# %%
brush.D(sigma=sigma, chi = chi_expected, N=300)
# %%
