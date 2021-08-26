from ascf_pb.solver import Pi, Phi, mu
def build_phi_profile_solver(kappa_cb, D_cb, phi_D_cb, **kwargs):
    kappa = kappa_cb(**kwargs)
    D = D_cb(kappa = kappa,**kwargs)
    phi_D = phi_D_cb(kappa = kappa, **kwargs)
    chi = kwargs['chi']
    def phi(z : float):
        return Phi(z, chi, kappa, D, phi_D)
    return phi, D

def build_Pi_profile_solver(*args, **kwargs):
    phi = build_phi_profile_solver(*args, **kwargs)[0]
    chi = kwargs['chi']
    def _Pi(z : float):
        return Pi(phi(z), chi)
    return _Pi

def get_D_mu(kappa_cb, D_cb, phi_D_cb, **kwargs):
    kappa = kappa_cb(**kwargs)
    D = D_cb(kappa = kappa,**kwargs)
    phi_D = phi_D_cb(kappa = kappa, **kwargs)
    chi = kwargs['chi']
    mu_D = mu(phi_D, chi)
    return D, mu_D

