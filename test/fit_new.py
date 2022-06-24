#%%
import ascf_pb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from scipy import optimize
#plt.rcParams["figure.dpi"] = 300
#plt.rcParams['figure.figsize'] = [12.0, 8.0]
LAST_USED_COLOR = lambda: plt.gca().lines[-1].get_color()


#Experimental data
Gamma = 4.9 #pmol/cm2
D_experimental_nm = 25 #nm
D_err_nm = 2 #nm
a = 0.76 #nm
N = 300
NA = 6.02214 #1e-23 mol-1
M=64.1 #kg/mol
conc = Gamma*M/D_experimental_nm*10 #mg/ml
rho=1400 #mg/ml
phi_experimental = conc/rho
#phi_err = N*sigma/(D_experimental)**2*D_err
phi_experimental_err = phi_experimental*D_err_nm/D_experimental_nm

#Expressed in SCF theory units
sigma = Gamma*a**2*NA*1e-3
D_experimental = D_experimental_nm/a
D_err = D_err_nm/a
D_min = N*sigma
#%%
brush = ascf_pb.BrushSolver()

#calculate chi on brush thickness
def find_chi_on_D(D, D_error, **kwargs):
    f = [
        lambda x: brush.D(chi = x, **kwargs)-D,
        lambda x: brush.D(chi = x, **kwargs)-D-D_error,
        lambda x: brush.D(chi = x, **kwargs)-D+D_error,
    ]
    chi = [optimize.brentq(f_, a=0.0, b=1.0) for f_ in f]
    return chi


chi = np.linspace(0, 1, 100)
D_on_chi = brush.D(N=N, sigma = sigma, chi = chi)
phi_av_on_chi = N*sigma/D_on_chi

#numerical derivative
d_chi = np.gradient(chi, D_on_chi)*a
chi_err = d_chi*D_err

chi_expected, chi_down, chi_up = find_chi_on_D(D_experimental, D_err, N=N, sigma=sigma)
phi_av = N*sigma/(D_experimental)
phi_err = N*sigma/(D_experimental)**2*D_err

#%%
annotation_text = f"D = {D_experimental_nm:.2f} $\pm$ {D_err_nm:.2f} nm" +\
    f"\n$\Gamma = {Gamma}$ pmol" +\
    f"\n$a = {a}$ nm" +\
    f"\n$N = {N}$" +\
    f"\n$\sigma = {sigma :.4f}$" +\
    "\n$\overline{\phi}_{theory} =$"+f"{phi_av :.3f}" + f"$\pm {phi_err :.3f}$\n" +\
    r"$\phi_{exp} = C/\rho =$"+ f"{phi_experimental :.3f}" + f"$\pm {phi_experimental_err :.3f}$\n" +\
    "$C =$"+f"{conc :.3f}"+ "mg/ml\n"
# %%
%matplotlib ipympl
fig, ax = plt.subplots()
ax.plot(D_on_chi*a, chi)
ax.plot(D_on_chi*a, chi-chi_err, linestyle = ":", color = LAST_USED_COLOR())
ax.plot(D_on_chi*a, chi+chi_err, linestyle = ":", color = LAST_USED_COLOR())

ax.set_ylabel("$\chi$")
ax.set_xlabel("D, [nm]")
ax.set_xlim(0,max(D_on_chi))

ax.axvline(x = D_min*a, linestyle = ":", color = "black")
ax.axhline(y = chi_expected, color = "red", linewidth = 0.2)
ax.scatter(x = D_experimental_nm, y=chi_expected, marker = "o", color = "red")
ax.vlines(x = D_experimental_nm, color = "red", ymin = chi_down, ymax = chi_up)
ax.text(x = D_experimental_nm, color = "red", y = chi_down-0.05, s = '{:.3f}'.format(chi_down), va = "top")
ax.text(x = D_experimental_nm, color = "red", y = chi_up+0.05, s = '{:.3f}'.format(chi_up), va = "bottom")

trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
ax.text(s="$D_{min} \equiv N\sigma / a$", x = D_min, y=0, va = "bottom", ha = "left", transform = trans)

ax.text(s=annotation_text, x = 0.99, y=0.99, va = "top", ha = "right", transform = ax.transAxes)

trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
ax.text(x = 0.98, color = "red", y = chi_expected, s = "$\chi =$"+'{:.3f}'.format(chi_expected), va = "bottom", ha = 'right', transform = trans)
# %%
fig, ax = plt.subplots()
ax.set_ylabel("$\chi$")
ax.set_xlabel("$\phi$")
ax.set_xlim(0,max(phi_av_on_chi))

ax.plot(phi_av_on_chi, chi)
ax.plot(phi_av_on_chi, chi-chi_err, linestyle = ":", color = LAST_USED_COLOR())
ax.plot(phi_av_on_chi, chi+chi_err, linestyle = ":", color = LAST_USED_COLOR())

ax.axhline(y = chi_expected, color = "red", linewidth = 0.2)
ax.scatter(x = phi_av, y=chi_expected, marker = "o", color = "red")
ax.vlines(x = phi_av, color = "red", ymin = chi_down, ymax = chi_up)
ax.text(x = phi_av, color = "red", y = chi_down-0.05, s = '{:.3f}'.format(chi_down), va = "top")
ax.text(x = phi_av, color = "red", y = chi_up+0.05, s = '{:.3f}'.format(chi_up), va = "bottom")

trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
ax.text(x = 0.98, color = "red", y = chi_expected, s = "$\chi =$"+'{:.3f}'.format(chi_expected), va = "bottom", ha = 'right', transform = trans)

ax.text(s=annotation_text, x = 0.01, y=0.99, va = "top", ha = "left", transform = ax.transAxes)
# %%
