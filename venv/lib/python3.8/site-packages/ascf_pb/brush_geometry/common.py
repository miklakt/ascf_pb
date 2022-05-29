from scipy.optimize import brentq
from ascf_pb.solver import Pi
from ascf_pb.solver import SolverError

def normalization_find_root(normalization, a, b, max_tries=20) -> float:
    tries = 0
    d=(b-a)/(max_tries+1)
    while tries < max_tries:
        try:
            root = brentq(normalization, a, b)
            break
        except SolverError as e:
            b = b-d
            print(f"brentq failed while normalization, changing bracketing interval")
            print(f"a:{a}, b:{b}")
            tries = tries + 1
    else:
        raise ValueError("Can't finnish normalization in several tries")
    return root

def phi_D_unrestricted(chi: float, **_) -> float:
    """Calculates polymer volume fraction at the end of the brush`
    Args:
        chi (float): Flory-Huggins parameter polymer-solvent
    Returns:
        float: volume fraction
    """
    almost_zero = 1e-06
    almost_one = 1.0 - almost_zero
    if chi <= 0.5:  # good solvent
        phi_D = 0
    else:  # poor solvent
        # to get phi_H we have to find where osmotic pressure vanishes
        def fsol(_phi):
            return Pi(_phi, chi)
        # find the root with brentq method
        try:
            phi_D = brentq(fsol, almost_zero, almost_one)
        except Exception as e:
            raise e
    return phi_D