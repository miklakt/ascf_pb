import numpy as np

from .solver import phi, Pi_phi

def mu_phi_chi(phi : float, chi : float) -> float:
    """Chemical potential for a given volume fraction and solvent regime

    Args:
        phi (float): local polymer volume fraction
        chi (float): Flory-Huggins parameter polymer-solvent

    Returns:
        float: chemical potential
    """    

    chem_pot = -np.log(1-phi)-2*chi*phi
    return chem_pot

def mu(sigma: float, chi: float, N: float):
    """Calculates chemical potential in a polymer brush at a given distance
    from the grafting surface

    Args:
        sigma (float): grafting density (chains per square)
        N (int): polymer chain's length
        chi (float): Flory-Huggins parameter polymer-solvent
        z (float): distance from the grafting surface
    Returns:
        function: chemical potential mu(z) function
    """
    _phi = phi(sigma, chi, N)
    @np.vectorize
    def _mu(z : float) -> float:
        """chemical potential mu(z) function compatible with numpy.array input
        Args:
            z (float): distance from the grafting surface
        Returns:
            float: chemical potential
        """
        return mu_phi_chi(_phi(z), chi)
    return _mu

def Pi(sigma: float, chi: float, N: float):
    """Calculates osmotic pressure in a polymer brush at a given distance
    from the grafting surface

    Args:
        sigma (float): grafting density (chains per square)
        N (int): polymer chain's length
        chi (float): Flory-Huggins parameter polymer-solvent
        z (float): distance from the grafting surface
    Returns:
        function: osmotic pressure Pi(z) function
    """
    _phi = phi(sigma, chi, N)
    @np.vectorize
    def _Pi(z : float) -> float:
        """osmotic pressure Pi(z) function for given sigma, N, chi
        function is compatible with numpy.array input
        Args:
            z (float): distance from the grafting surface
        Returns:
            float: osmotic pressure
        """
        return Pi_phi(_phi(z), chi)
    return _Pi