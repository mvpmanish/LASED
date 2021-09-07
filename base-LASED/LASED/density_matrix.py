'''
Density matrix operations
Author: Manish Patel
Date Created: 11/06/2021
'''

import numpy as np
from LASED.index import *

def getSingleStateMatrix(flat_rho, n, sub_states):
    """Obtain an angular momentum state density matrix from the flattened coupled state density rho vector.
    
    Parameters:
    flat_rho : array of arrays with one column of all density matrix elements of coupled E & G states
    n : number of states in total laser-coupled system
    sub_states : a list of the excited or ground states, E or G respectively
    
    Returns:
    A square matrix of size length of sub_states. Elements are ordered from left to right
    according to the order of sub_states e.g. if state labelled 1 is first then first element
    would correspond to rho_11 i.e. population of state 1
    """
    # Set up the state matrix
    state_density_matrix = np.zeros((len(sub_states), len(sub_states)), dtype = np.complex)
    
    # Populate matrix
    for i, sub_state in enumerate(sub_states):
        for j, sub_state_p in enumerate(sub_states):
            state_density_matrix[i, j] = flat_rho[index(sub_state, sub_state_p, n), 0]
    return state_density_matrix


def JNumber(state_list):
    """
    Calculate the angular momentum from the number of states in the list.
    """
    return int((len(state_list)-1)/2)  # J = (m-1)/2


def appendDensityMatrixToFlatCoupledMatrix(flatrho, density_rho, sub_states, n):
    """ Adds density matrix elements to a flat, coupled matrix.
    
    Parameters:
    - flat_rho: array of arrays with one column of all density matrix elements of coupled E and G states
    - n: number of states in total laser-coupled system
    - sub-states: a list of the excited or ground states
    - density_rho: either the excited or ground state density matrix with the convention that the 
                   upper left-hand of the matrix is state population for m_J = -J if the state has angular momentum J 
"""
    for i, sub_state in enumerate(sub_states):
        for j, sub_state_p in enumerate(sub_states):
            flatrho[index(sub_state, sub_state_p, n), 0] = density_rho[i, j]         
