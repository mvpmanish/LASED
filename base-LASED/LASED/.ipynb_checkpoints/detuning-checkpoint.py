'''
Define functions for the detuning of an atomic system
'''
from LASED.state import *
from LASED.constants import *

def delta(e, g):
    """Detunings between substates.
     
     Parameters:
        e (State): State object
        g (State): State object
    
    Returns:
        float: Difference in angular frequency of states (Grad/s).
    """
    return e.w - g.w

def angularFreq(wavelength):
    """Calculates the angular frequency in Grad/s from a given wavelength.
    
    Parameters:
        wavelength (float): A wavelength in nm
    
    Returns:
        float: The angular frequency in Grad/s
    """
    return 2*PI*C/wavelength*1e-9
    
def dopplerDelta(e, g, w_q, lambda_q, v_z):
    """ The detuning between excited and ground states.
    
    Accounts for a fixed motion of the atoms. Used between excited and ground states.
    
    Parameters:
        e (State): State object for excited state.
        g (State): State object for ground state.
        w_q (float): Angular frequency of exciting laser in rad/s.
        lambda_q (float): Wavelength of exciting laser in m.
        v_z (float): Velocity component of atoms in direction of laser in m/s.
    
    Returns:
        float: The detuning between ground and excited states including the doppler detuning due to a given atomic velocity.
    
    """
    return w_q - v_z/lambda_q - e.w + g.w