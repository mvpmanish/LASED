'''
Define functions for the detuning of an atomic system
'''
from state import *
from constants import *

# Detunings
# Input: e is a State object, g is a State object
def delta(e, g):
    """
    Detunings between states. 
    Inputs:
        e: State object
        g: State object
    Output:
        (float) difference in angular frequency of states (Grad/s)
    """
    return e.w - g.w

def angularFreq(wavelength):
    """
    Calculates the angular frequency in Grad/s from a given wavelength.
    Inputs:
        wavelength (float): a wavelength in nm
    Output:
        (float) the angular frequency in Grad/s
    """
    return 2*PI*C/wavelength*1e-9
    
# Doppler detuning.
# Input: w_q is the angular frequency of the laser in rad/s, lambda_q is the wavelength of light in m
# ,v_z is the velocity component of atoms in direction of laser in m/s, e and g are State objects
def dopplerDelta(e, g, w_q, lambda_q, v_z):
    return w_q - v_z/lambda_q - e.w + g.w