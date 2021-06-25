'''
Definition of the half-Rabi frequency and all other functions needed to calculate it
'''

import numpy as np
from state import *
from constants import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j

# Intensity of Gaussian TEM_00 laser beam profile
# Inputs: P_las is power of laser given by power meter in mW, r_sigma is radius at 2D std in mm, r is radius in mm
# Output: intensity at radius r from beam axis given a Gaussian intesnity profile
def gaussianIntensity(P_las, r_sigma, r):
    return (P_las/(2*PI*r_sigma*r_sigma))*np.exp(-(r*r)/(2*r_sigma*r_sigma))

# Polarisation quantum number
def quantumPolarisation(q):
    if q == 1:
        return 0
    if (q == -1) or (q == 0):
        return 1

# Coupling constants
# Input: e and g are State objects, q is the laser polarisation
# Output: the coupling coefficient which is the Dipole matrix element
def coupling(e, g, q):
    sign = np.power(-1, quantumPolarisation(q)+e.F+g.F+e.J+g.J+e.I+e.L+e.S-e.m+ 1)
    factor = np.sqrt((2*e.F+1)*(2*g.F+1)*(2*e.J+1)*(2*g.J+1)*(2*e.L+1))
    wig6j_1 = wigner_6j(e.J, e.F, e.I, g.J, g.L, 1)
    wig6j_2 = wigner_6j(e.L, e.J, e.S, g.J, g.L, 1)
    wig3j = wigner_3j(e.F, 1, g.F, -1*e.m, q, g.m)
    return sign*factor*wig6j_1*wig6j_2*wig3j

# Caluclates the half-Rabi frequency in GHz.
# Inputs: intensity given in mW/mm^2, lifetime in ns/rad, wavelength of laser in m
# Output: the half-Rabi frequency constant (without coupling coefficient) in GHz
def halfRabiFreq(intensity, lifetime, wavelength):
    I = intensity*1000  # Convert mW/mm^2 to W/m^2
    h = 6.6260715e-34  # Plank's constant
    c = 299792458  # Speed of light
    tau = lifetime*1e-9  # Convert ns/rad lfetime to s/rad
    
    return np.sqrt((3*I*wavelength**3)/(8*np.pi*c*h*tau))*1e-9  #Gives half-Rabi freq in Grad/s