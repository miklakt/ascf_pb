from .profiles import build_phi_profile_solver,  build_Pi_profile_solver
from .topology.plain import D_unrestricted, phi_D_unrestricted
from .topology.kappa import kappa_plain
from functools import partial
from numpy import vectorize


def D(**kwargs):
    kappa = kappa_plain(kwargs['N'])
    return D_unrestricted(kappa=kappa, **kwargs)

phi = partial(build_phi_profile_solver,
    kappa_cb=kappa_plain, D_cb = D_unrestricted, 
    phi_D_cb=phi_D_unrestricted)

Pi = partial(build_phi_profile_solver,
    kappa_cb=kappa_plain, D_cb = D_unrestricted, 
    phi_D_cb=phi_D_unrestricted)

if __name__ == '__main__':
    import numpy as np
    args = dict(chi=0, N=1000, sigma = 0.02)
    z = np.arange(round(D(**args)+0.5))
    print([phi(**args)(z_) for z_ in z])