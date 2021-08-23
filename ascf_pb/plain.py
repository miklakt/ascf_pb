from solver import Pi, Phi
from scipy.optimize import brentq
from scipy import integrate
import numpy as np

def phi_D_unrestricted(chi: float) -> float:
    """Calculates polymer volume fraction at the end of the brush`
    Args:
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: volume fraction
    """
    if chi <= 0.5:  # good solvent
        phi_D = 0
    else:  # poor solvent
        # to get phi_H we have to find where osmotic pressure vanishes
        def fsol(_phi):
            return Pi(_phi, chi)
        # find the root with brentq method
        min_phi = 0.00001  # exclude 0 from the roots
        max_phi = 0.99999  # exclude 1 from the roots
        try:
            phi_D = brentq(fsol, min_phi, max_phi)
        except Exception as e:
            raise e
    return phi_D

def normalization_plain_unrestricted(chi : float, kappa : float, theta : float):
    phi_D  = phi_D_unrestricted(chi)
    def integrand(z, d):
        return Phi(z = z, d=d,
        chi = chi, kappa=kappa, phi_D=phi_D)
    def integral(d):
        return integrate.quad(integrand, 0, d, args=(d,))[0] - theta
    return integral

def D_unrestricted(chi : float, kappa : float, theta : float, max_D_guess : float):
    normalization = normalization_plain_unrestricted(chi,kappa,theta)
    min_D = 0
    max_tries = 20
    tries = 0
    dD=max_D_guess/(max_tries+1)
    while tries < max_tries:
        try:

            _D = brentq(normalization, min_D, max_D_guess)
            break
        except ValueError as e:
            #print(e)
            #print(f'trying to decrease upper boundary max_H_guess')
            max_D_guess = max_D_guess-dD
            #print(f"max_H_guess : {max_H_guess}")
            tries = tries + 1
    else:
        raise ValueError()
    return _D

def D(chi, N, sigma, kappa = None):
    if kappa is None:
        from ascf_pb.kappa import kappa_plain
        kappa = kappa_plain(N)
    theta = N*sigma
    return D_unrestricted(chi, kappa, theta, N)

def phi(chi, N, sigma, kappa = None):
    if kappa is None:
        from ascf_pb.kappa import kappa_plain
        kappa = kappa_plain(N)
    _D = D(chi, N, sigma, kappa)
    phi_D = phi_D_unrestricted(chi)
    @np.vectorize
    def phi(z : float):
        return Phi(z, chi, kappa, _D, phi_D)
    return phi

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""
    Planar polymer brush profiles calculation 
    using Analytic Self-Consistent Field method
    """
                                        )
    parser.add_argument('sigma', type = float,
                        help = 'grafting density (chains per square)')
    parser.add_argument('chi', type = float,
                        help = 'Flory-Huggins parameter polymer-solvent')
    parser.add_argument('N', type = int,
                        help = "polymer chain's length")

    args = parser.parse_args()

    _D = D(args.chi, args.N, args.sigma)
    print('D:', _D)
    z  = np.arange(round(_D+1.5))
    _phi = phi(args.chi, args.N, args.sigma)(z)
    print('phi:', _phi)

