#%%
from sympy import solve, Symbol, Eq, series, log, sqrt
chi = Symbol('chi', positive=True)
phi= Symbol('phi', positive=True)
mu_known = Symbol('mu')
# %%
mu = -log(1-phi)+2*chi*phi
# %%
phi_solve = solve(Eq(mu_known, mu),phi)[0]
# %%
phi_solve
#%%
Lambda = Symbol('Lambda')
z = Symbol('z')
k = Symbol('k')
H = Symbol('H')
phi_D = Symbol('phi_D')
mu_D = -log(1-phi_D)+2*chi*phi_solve

z = Symbol('z')
#%%
phi_z = phi_solve.subs(mu_known, k*(Lambda**2-z**2))
# %%
from sympy import integrate, exp
norm = integrate(phi_solve, mu_known)

# %%
norm
# %%
