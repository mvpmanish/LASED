'''
Definition of the half-Rabi frequency and all other functions needed to calculate it
'''

import numpy as np
from LASED.state import *
from LASED.constants import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j

# 
# Inputs: P_las is , r_sigma is , r is 
# Output: 
def gaussianIntensity(P_las, r_sigma, r):
    """Intensity of Gaussian TEM_00 laser beam profile.
    
    Parameters:
        P_las (float): Power of laser given by power meter in mW
        r_sigma (float): Radius at 2D standard deviation of a Gaussian in mm
        r (float): Radius in mm at which the laser intensity is evaluated
        
    Returns:
        Intensity in mW/mm^2 at radius r from beam axis given a Gaussian intensity profile
    """
    return (P_las/(2*PI*r_sigma*r_sigma))*np.exp(-(r*r)/(2*r_sigma*r_sigma))

# 
# Input: e and g are State objects, q is the laser polarisation
# Output: the coupling coefficient which is the Dipole matrix element
def coupling(e, g, q):
    """Coupling constants between states via laser radiation
    
    Parameters: 
        e (State): the excited state
        g (State): the ground state
        q (int): the laser polarisation coupling e and g. This can be -1 for LH circular, +1 for RH circular, or 0 for linear polarisation.
        
        Returns:
            (float) the coupling constant between an excited and ground state
    """
    sign = np.power(-1, (q*(1+q)/2)+e.F+g.F+e.J+g.J+e.I+e.L+e.S-e.m+1)
    factor = np.sqrt((2*e.F+1)*(2*g.F+1)*(2*e.J+1)*(2*g.J+1)*(2*e.L+1))
    wig6j_1 = wigner_6j(e.J, e.F, e.I, g.F, g.J, 1)
    wig6j_2 = wigner_6j(e.L, e.J, e.S, g.J, g.L, 1)
    wig3j = wigner_3j(e.F, 1, g.F, -1*e.m, q, g.m)
    return sign*factor*wig6j_1*wig6j_2*wig3j

# 
# Inputs: intensity given in mW/mm^2, lifetime in ns, wavelength of laser in m
# Output: the half-Rabi frequency constant (without coupling coefficient) in GHz
def halfRabiFreq(intensity, lifetime, wavelength):
    """Calculates the half-Rabi frequency in Grad/s.
     
    Parameters: 
        intensity: intensity of laser in mW/mm^2
        lifetime: lifetime of the excited state to the ground state transition in nanoseconds
        wavelength: wavelength (in metres) corresponding to the resonant transition from the ground to excited state.
    
    Returns:
        (float) the half-Rabi frequency in Grad/s
    
    """
    I = intensity*1000  # Convert mW/mm^2 to W/m^2
    h = 6.6260715e-34  # Plank's constant
    #hbar = 1.0456e-34
    c = 299792458  # Speed of light
    tau = lifetime*1e-9  # Convert ns lfetime to s
    
    return np.sqrt((3*I*wavelength**3)/(8*PI*c*h*tau))*1e-9  #Gives half-Rabi freq in Grad/s