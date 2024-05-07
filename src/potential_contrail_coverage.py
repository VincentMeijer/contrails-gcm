import numpy as np

def get_cirrus_coverage(RHi, RHi_ci=0.6, RHi_sat=1.0):
    """
    Implements Equation 1 in Burkhardt et al. (2008).

    Parameters
    ----------
    RHi : float
        Relative humidity w.r.t. ice
    RHi_ci : float
        Relative humidity w.r.t. ice above which part of the grid box will 
        form a cirrus cloud
    RHi_sat : float
        Relative humidity w.r.t. ice at which the full grid box is cloudy

    Returns
    -------
    b_ci : float
        The fraction of the grid box that is cloudy
    """
    
    return 1 - (1 - (RHi - RHi_ci) / (RHi_sat - RHi_ci))**0.5


def get_RHi_nuc(T):
    """
    Returns the homogeneous freezing threshold for cirrus according to
    Koop (2004).

    Parameters
    ----------
    T : float
        Temperature in Kelvin

    Returns
    -------
    RHi_nuc : float
        Homogeneous freezing threshold for cirrus
    """
    return 2.349 - T / 259

def get_RHi_star(RHi_sat, RHi_ci, RHi_cc):
    """
    Parameters
    ----------
    RHi_sat : float
        Relative humidity w.r.t. ice at which the full grid box is cloudy
    RHi_ci : float
        Relative humidity w.r.t. ice above which part of the grid box will 
        form a cirrus cloud
    RHi_cc : float
        Relative humidity w.r.t. ice above which part of the grid box will
        be ice supersaturated

    Returns
    -------
    RHi_star : float
    """
    return RHi_sat - (RHi_ci - RHi_cc)**2 / (RHi_sat - RHi_ci)

def get_contrail_cirrus_coverage(RHi, RHi_ci=0.6, RHi_cc=0.4, RHi_sat=1.0):
    
    b_ci = np.maximum(0, get_cirrus_coverage(RHi, RHi_ci=RHi_ci,
                                             RHi_sat=RHi_sat))
    RHi_star = get_RHi_star(RHi_sat, RHi_ci, RHi_cc)

    B_cc_ci = (RHi - RHi_cc) / (RHi_sat - RHi_ci) - b_ci * (1 - b_ci)
    B_cc_ci[RHi > RHi_star] = 1.0
    return B_cc_ci
    

    