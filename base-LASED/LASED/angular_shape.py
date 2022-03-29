'''Define functions for creating an angular shape of a single state for an atom and the time evolution of the shape.
Author: Manish Patel
Date created: 24/03/2022
'''

from LASED.density_matrix import *
import numpy as np
from scipy.special import sph_harm

def getFlattenedRhot(rho_t, n):
    """Flatten a density matrix over time rho(t).
    
    Parameters:
        rho_t (list of list): List of lists of time evolution of density matrix elements.
        n (int): Number of states in total laser-coupled system.
        
    Returns:
        (List of (list of lists)): List of flattened density matrices for the time evolution.
    """
    flat_rho_t = []
    for rho in np.transpose(rho_t):
        new_rho = np.zeros((n*n, 1), dtype = complex)  # Placeholder
        for i, element in enumerate(new_rho):
            new_rho[i, 0] = rho[i]
        flat_rho_t.append(new_rho)
    return flat_rho_t

def angularShape(flat_rho, n, sub_states, theta, phi):
    """Computes the angular shape of the single atomic state |J, m>.
    
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
    J = sub_states[0].F  # Get the total angular momentum
    rho = getSingleStateMatrix(flat_rho, n, sub_states)  # Get the single state density matrix
    W = []
    for thta in theta:
        Wrow = []
        for ph in phi:
            Wsum = 0
            for i, sub_state in enumerate(sub_states):
                for j, sub_state_p in enumerate(sub_states):
                    Wsum += (rho[i, j]*sph_harm(sub_state.m, J, thta, ph)*sph_harm(sub_state_p.m, J, thta, ph).conjugate())
            Wrow.append(abs(Wsum))
        W.append(Wrow)
    return W

def SphericalToCartesian(r, theta, phi):
    """Turns spherical coordinates to cartesian coordianates.
    
    Parameters:
        r (float): Radius from 0 to infinity.
        theta (array_like): Azimuthal (longitudinal) coordinate in [0, 2*pi].
        phi (array_like): Polar (colatitudinal) coordinate in [0, pi].
        
    Returns:
        x,y,z (tuple): The cartesian coordinates.
    """
    x = r* np.cos(theta)* np.sin(phi)
    y = r* np.sin(theta)* np.sin(phi)
    z = r* np.cos(phi)
    return x, y, z
