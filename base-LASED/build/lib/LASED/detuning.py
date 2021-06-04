'''
Define functions for the detuning of an atomic system
'''
from state import *

# Detunings
# Input: e is a State object, g is a State object
def Delta(e, g):
    return e.w - g.w
    
# Doppler detuning.
# Input: w_q is the angular frequency of the laser in rad/s, lambda_q is the wavelength of light in m
# ,v_z is the velocity component of atoms in direction of laser, e and g are State objects
def dopplerDelta(e, g, w_q, lambda_q, v_z):
    return w_q - v_z/lambda_q + e.w - g.w