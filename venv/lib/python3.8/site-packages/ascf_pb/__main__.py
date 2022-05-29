if __name__=="__main__":
    import argparse
    import ascf_pb
    import sys
    import os
    import contextlib
    import numpy as np
    parser = argparse.ArgumentParser(description='IP')
    parser.add_argument('N', 
                        metavar='N',
                        type = float,
                        help = 'the chain length of a polymer brush',
                        #required=True
                        )
    parser.add_argument('sigma', 
                        metavar='sigma',
                        type = float,
                        help = 'number of chains per unit area',
                        #required=True
                        )
    parser.add_argument('chi', 
                        metavar='chi',
                        type = float,
                        help = 'Flory-Huggins polymer-solvent interaction parameter',
                        #required=True
                        )
    
    parser.add_argument('-eta',
                        metavar='eta', 
                        help = 'topological parameter default is 1.0',
                        type = float,
                        required=False,
                        default = 1.0)

    parser.add_argument('-R',
                        metavar='R',
                        help = 'distance to the restriction surface',
                        type = float,
                        required=False)
    
    parser.add_argument('-pore_R', 
                        metavar='pore_Radius',
                        help = 'pore radius',
                        type = float,
                        required=False)

    parser.add_argument('-t', 
                        metavar='topology',
                        help = 'system geometry: plain, pore',
                        type = str,
                        required=False,
                        default = 'plain')

    parser.add_argument('-z_step', 
                        metavar='z_step',
                        help = 'profile fidelity',
                        type = float,
                        required=False,
                        default = 1.0)

    args = parser.parse_args()


    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
        D = ascf_pb.D(**vars(args))()
        z = np.arange(0, int(D/args.z_step+1.5)*args.z_step, args.z_step)
        phi_func = np.vectorize(ascf_pb.phi(**vars(args)))
        phi = phi_func(z)
    print(D)
    print(z)
    print(phi)
    
