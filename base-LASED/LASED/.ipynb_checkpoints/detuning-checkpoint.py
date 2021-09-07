'''
Define functions for the detuning of an atomic system
'''
from LASED.state import *
from LASED.constants import *

def delta(e, g):
    """Detunings between substates.
     
     Parameters:
        e: State object
        g: State object
    
    Returns:
        (float) difference in angular frequency of states (Grad/s)
    """
    return e.w - g.w

def angularFreq(wavelength):
    """Calculates the angular frequency in Grad/s from a given wavelength.
    
    Parameters:
        wavelength (float): a wavelength in nm
    
    Returns:
        (float) the angular frequency in Grad/s
    """
    return 2*PI*C/wavelength*1e-9
    
def dopplerDelta(e, g, w_q, lambda_q, v_z):
    """ The detuning between excited and ground states.
    
    Accounts for a fixed motion of the atoms. Used between excited and ground states.
    
    Parameters:
        e: State object for excited state
        g: State object for ground state
        w_q: Angular frequency of exciting laser in rad/s
        lambda_q: Wavelength of exciting laser in m
        v_z: Velocity component of atoms in direction of laser in m/s
    """
    return w_q - v_z/lambda_q - e.w + g.w