from .solver import Pi, Phi
def build_phi_profile_solver(kappa_cb, D_cb, phi_D_cb, **kwargs):
    kappa = kappa_cb(**kwargs)
    D = D_cb(kappa = kappa,**kwargs)
    phi_D = phi_D_cb(**kwargs)
    chi = kwargs['chi']
    def phi(z : float):
        return Phi(z, chi, kappa, D, phi_D)
    return phi

def build_Pi_profile_solver(*args, **kwargs):
    phi = build_phi_profile_solver(*args, **kwargs)
    chi = kwargs['chi']
    def _Pi(z : float):
        return phi(z)
    return _Pi

