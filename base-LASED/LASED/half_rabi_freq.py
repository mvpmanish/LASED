'''
Definition of the half-Rabi frequency and all other functions needed to calculate it
'''

import numpy as np
from LASED.state import *
from LASED.constants import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from sympy import Rational
from functools import lru_cache

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
        float: Intensity in mW/mm^2 at radius r from beam axis given a Gaussian intensity profile
    """
    return (P_las/(2*PI*r_sigma*r_sigma))*np.exp(-(r*r)/(2*r_sigma*r_sigma))

def _is_triangle(a, b, c):
    return abs(a - b) <= c <= a + b


@lru_cache(maxsize=None)
def _coupling_cached(eL, eS, eJ, eI, eF, em, gL, gS, gJ, gF, gm, q):
    # Arguments are twice the quantum numbers (exact integers); q is the polarisation.
    # Selection rules: skip the (slow, exact) sympy evaluation when the result is zero
    if -em + 2*q + gm != 0 or not _is_triangle(eF, gF, 2) or not _is_triangle(eJ, gJ, 2) \
            or not _is_triangle(eL, gL, 2):
        return 0.0
    eL, eS, eJ, eI, eF, em, gL, gS, gJ, gF, gm = (Rational(x, 2) for x in (eL, eS, eJ, eI, eF, em, gL, gS, gJ, gF, gm))
    exponent = (q*(1+q)/2)+eF+gF+eJ+gJ+eI+eL+eS-em+1
    sign = -1.0 if int(round(exponent)) % 2 else 1.0
    factor = np.sqrt(float((2*eF+1)*(2*gF+1)*(2*eJ+1)*(2*gJ+1)*(2*eL+1)))
    wig6j_1 = wigner_6j(eJ, eF, eI, gF, gJ, 1)
    wig6j_2 = wigner_6j(eL, eJ, eS, gJ, gL, 1)
    wig3j = wigner_3j(eF, 1, gF, -1*em, q, gm)
    return float(sign*factor*wig6j_1*wig6j_2*wig3j)


def coupling(e, g, q):
    """Coupling constants between states via laser radiation

    This is the (dimensionless) dipole matrix element <e| d.u_q |g> between the sub-states in units of the
    reduced matrix element, with u_{+1} = (x + iy)/sqrt(2), u_0 = z, u_{-1} = (x - iy)/sqrt(2).
    Results are cached, so repeated calls are cheap.

    Parameters:
        e (State): the excited state
        g (State): the ground state
        q (int): the laser polarisation coupling e and g. This can be -1 for LH circular, +1 for RH circular, or 0 for linear polarisation.

        Returns:
            float: the coupling constant between an excited and ground state
    """
    return _coupling_cached(_twice(e.L), _twice(e.S), _twice(e.J), _twice(e.I), _twice(e.F), _twice(e.m),
                            _twice(g.L), _twice(g.S), _twice(g.J), _twice(g.F), _twice(g.m), int(q))


def _twice(x):
    """Twice an integer or half-integer quantum number (possibly given as a float), as an exact int."""
    r = int(round(2*float(x)))
    if abs(r - 2*float(x)) > 1e-9:
        raise ValueError(f"Quantum number {x} is not an integer or half-integer.")
    return r


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
        float: The half-Rabi frequency in Grad/s.
    
    """
    I = intensity*1000  # Convert mW/mm^2 to W/m^2
    h = 6.6260715e-34  # Plank's constant
    #hbar = 1.0456e-34
    c = 299792458  # Speed of light
    tau = lifetime*1e-9  # Convert ns lfetime to s
    
    return np.sqrt((3*I*wavelength**3)/(8*PI*c*h*tau))*1e-9  #Gives half-Rabi freq in Grad/s