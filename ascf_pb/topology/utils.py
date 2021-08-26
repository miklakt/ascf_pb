from scipy.optimize import brentq
def normalization_find_root(normalization, a, b, max_tries=20):
    tries = 0
    d=b/(max_tries+1)
    while tries < max_tries:
        try:
            root = brentq(normalization, a, b)
            break
        except ValueError as e:
            b = b-d
            print(f"brentq failed while normalization, changing bracketing interval")
            tries = tries + 1
    else:
        raise ValueError("f(a) and f(b) still have the same signs")
    return root