# __main.py__

import argparse
import numpy as np

import ascf_pb

def test():
    _H = ascf_pb.H(sigma=0.05, chi=0, N=1000)
    print(_H)

def main():
    
    def called_from_shell(sigma: float, chi: float, N: int):
        _H = ascf_pb.H(sigma, chi, N)
        print(f"\nH : {_H}")
        z = np.arange(0, np.ceil(_H))
        _phi = ascf_pb.phi(sigma,chi,N)(z)
        _Pi = ascf_pb.Pi(sigma,chi,N)(z)
        _mu = ascf_pb.mu(sigma,chi,N)(z)
        print("\nphi :")
        print(*_phi, sep = '\n')
        print("\nPi :")
        print(*_Pi, sep = '\n')
        print("\nmu :")
        print(*_mu, sep = '\n')


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

    args = vars(parser.parse_args())

    called_from_shell(**args)

if __name__ == "__main__":
    test()