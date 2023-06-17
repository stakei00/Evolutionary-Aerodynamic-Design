import numpy as np 
import scipy.interpolate

"""
below are some preset fitness functions which can be selected from
"""

def max_LtoD(wing:object) -> float:
    """
    fitness based on the highest lift-to-drag ratio
    This corresponds to the best range CL for a propeller-driven aircraft and 
    best endurance CL for a jet-powered aircraft 
    """
    CL_interp = np.linspace(min(wing.lift_coefficient), max(wing.lift_coefficient), 100)
    CD_interp = np.interp(CL_interp, wing.lift_coefficient, wing.drag_coefficient)
    lift_to_drag = np.divide(CL_interp, CD_interp)
    return max(lift_to_drag)


def max_L12toD(wing:object) -> float: 
    """
    fitness based on the highest sqrt(lift)/drag ratio
    This corresponds to the best range CL for a jet-powered aircraft 
    """
    CL_interp = np.linspace(min(wing.lift_coefficient), max(wing.lift_coefficient), 100)
    CD_interp = np.interp(CL_interp, wing.lift_coefficient, wing.drag_coefficient)
    lift12_to_drag = np.divide(np.power(CL_interp,1/2), CD_interp)
    return max(lift12_to_drag)


def LtoD_at_CL(wing:object, CL:float=None) -> float:
    """
    fitness based on the highest lift-to-drag ratio at specified lift coefficient.
    Use this if cruise/loiter CL is already known ahead of time.
    """
    lift_to_drag = np.divide(wing.lift_coefficient,wing.drag_coefficient)
    L2D_interp_func = scipy.interpolate.interp1d(wing.lift_coefficient, lift_to_drag, kind="cubic")
    return L2D_interp_func(CL)


def max_L32toD(wing:object) -> float:
    """
    fitness based on highest L^(3/2)/D ratio
    This corresponds the best endurace CL for a propeller-driven aircraft 
    """
    CL_interp = np.linspace(min(wing.lift_coefficient), max(wing.lift_coefficient), 100)
    CD_interp = np.interp(CL_interp, wing.lift_coefficient, wing.drag_coefficient)
    lift32_to_drag = np.divide(np.power(CL_interp,3/2), CD_interp)
    return max(lift32_to_drag)


def max_CL_max(wing:object) -> float:
    """
    fitness based on highest maximum lift coefficient. Use this for designing 
    very-very slow flying aircraft 
    """
    return wing.max_lift_coefficient 