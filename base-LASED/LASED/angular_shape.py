'''Define functions for creating an angular shape of a single state for an atom and the time evolution of the shape.
Author: Manish Patel
Date created: 24/03/2022
'''

from LASED.density_matrix import *
import numpy as np
from scipy.special import sph_harm

def angularShape(flat_rho, n, sub_states, J, theta, phi):
    """ Computes the angular shape of the single atomic state |J, m>.
    
    Parameters:
        flat_rho (list of list): Flattened 2D density matrix of coupled E & G states with the laser.
        n (int): Number of states in total laser-coupled system.
        sub_states (list): List of sub_states with angular momentum J.
        J (int): Total angular momentum quantum number
        theta (array_like): Azimuthal (longitudinal) coordinate in [0, 2*pi].
        phi (array_like): Polar (colatitudinal) coordinate in [0, pi].
    
    Returns:
        (array_like): The radius of the angular shape of the input atomic state |J, m>.
    """
    
    getSingleStateMatrix(flat_rho, n, sub_states)  # Get the single state density matrix
    W = []
    for ph in phi:
        Wrow = []
        for thta in theta:
            Wsum = 0
            for i, sub_state in enumerate(sub_states):
                for j, sub_state_p in enumerate(sub_states):
                    Wsum += (rho[i, j]*sph_harm(sub_state.m, J, thta, ph)*sph_harm(sub_state_p.m, J, thta, ph).conjugate())
            Wrow.append(abs(Wsum))
        W.append(Wrow)
    return W
