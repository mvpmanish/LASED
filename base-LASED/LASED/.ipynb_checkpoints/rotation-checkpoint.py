'''
This is a file to define functions for rotating density matrices in the QED simulation of
a laser-atom system.
Author: Manish Patel
Date created: 12/05/2021
'''

import numpy as np
import math

def wigner_D(J, alpha, beta, gamma):
    '''
    Calculates the Wigner D-matrix for rotation by Eueler angles (alpha, beta, gamma).
    Inputs:
        J: total angular momentum quantum number of the state which will be rotated with the 
        resulting D-matrix
        alpha: rotation around z-axis
        beta: rotation about the y'-axis
        gamma: rotation about the z''-axis
    Returns:
        A square matrix of size 2J+1
    '''
    size = 2*J+1  # Number of sub-states
    m = np.linspace(-J, J, size, dtype=int)  # Projections of J
    D = np.zeros((size, size), dtype = np.complex)  # Set up D-matrix
    for i, mp in enumerate(m):
        for j, mpp in enumerate(m):
            print("i, j:", i, j)
            print("mp, mpp:", mp, mpp)
            alpha_const = math.cos(-mp*alpha)+1.j*math.sin(-mp*alpha)
            gamma_const = math.cos(-mpp*gamma)+1.j*math.sin(-mpp*gamma)
            D[i, j] = alpha_const*small_Wigner_D(J, beta, mp, mpp)*gamma_const
    return D

def small_Wigner_D(J, beta, mp, m):
    '''
    Calculates the small Wigner D-matrix elements for rotation by Euler angles (alpha, beta, gamma)
    '''
    const = np.sqrt((math.factorial(J+mp))*math.factorial(J-mp)*math.factorial(J+m)*math.factorial(J-m))
    print(const)
    d_sum = 0
    # Define limits so sum does not contain negative factorials
    s_max = min(J+m, J-mp)
    s_min = max(0, m-mp)
    sum_index = np.linspace(s_min, s_max, )
    for s in range(s_min, s_max+1):  # Have to go to s_max+1 or will miss out on the s_max value
        numerator = np.power(-1, mp - m + s)*np.power(math.cos(beta/2), 2*J+m-mp-2*s)*np.power(math.sin(beta/2), mp-m+2*s)
        denominator = math.factorial(J+m-s)*math.factorial(s)*math.factorial(mp-m+s)*math.factorial(J-mp-s)
        d_sum += numerator/denominator
    
    return const*d_sum

'''
Rotates a density matrix by the wigner D-matrix determined by Euler angles (alpha, beta, gamma).
(rho_newframe) = (D)(rho)(D*). The z-y-z convention is used where alpha rotates around the z-axis,
then beta rotates around the new y-axis, and gamma rotates around the new z-axis
'''
def rotation(rho, J, alpha, beta, gamma):
    D_matrix = wigner_D(J, alpha, beta, gamma)
    D_conj = np.transpose(np.conj(D_matrix))
    return np.dot(D_matrix, np.dot(rho, D_conj))

'''
Obtain an angular momentum state matrix from the flattened coupled state density rho vector
Inputs:
    flat_rho : array of arrays with one column of all density matrix elements of coupled E & G states
    n : number of states in total laser-coupled system
    sub_states : a list of the excited or ground states, E or G respectively
Returns:
    A square matrix of size length of sub_states. Elements are ordered from left to right
    according to the order of sub_states e.g. if state labelled 1 is first then first element
    would correspond to rho_11 i.e. population of state 1
'''
def getSingleStateMatrix(flat_rho, n, sub_states):
    # Set up the state matrix
    state_density_matrix = np.zeros((len(sub_states), len(sub_states)), dtype = np.complex)
    
    # Populate matrix
    for i, sub_state in enumerate(sub_states):
        for j, sub_state_p in enumerate(sub_states):
            state_density_matrix[i, j] = flat_rho[index(sub_state, sub_state_p, n), 0]
    return state_density_matrix
'''
Calculate the angular momentum from the number of states in the list
'''
def JNumber(state_list):
    return int((len(state_list)-1)/2)  # J = (m-1)/2

'''
Adds density matrix elements to a flat, coupled matrix.
Inputs:
    - flat_rho: array of arrays with one column of all density matrix elements of coupled E and G states
    - n: number of states in total laser-coupled system
    - sub-states: a list of the excited or ground states
    - density_rho: either the excited or ground state density matrix with the convention that the 
                   upper left-hand of the matrix is state population for m_J = -J if the state has angular momentum J 
''' 
def appendDensityMatrixToFlatCoupledMatrix(flatrho, density_rho, sub_states, n):
    for i, sub_state in enumerate(sub_states):
        for j, sub_state_p in enumerate(sub_states):
            flatrho[index(sub_state, sub_state_p, n), 0] = density_rho[i, j]         

'''
Rotate the excited and ground state populations by the Euler angles alpha, beta, gamma
'''
def rotateInitialMatrix(flat_rho, n, E, G, alpha, beta, gamma):
    # Make a copy to return
    rotated_rho = copy.deepcopy(flat_rho)
    # Rotate the excited state populations and atomic coherences
    J_E = JNumber(E)
    rotated_excited_rho = rotation(getSingleStateMatrix(rotated_rho, n, E), J_E, alpha, beta, gamma)
    appendDensityMatrixToFlatCoupledMatrix(rotated_rho, rotated_excited_rho, E, n)
    
    # Rotate the ground state populations and atomic coherences
    J_G = JNumber(G)
    rotated_ground_rho = rotation(getSingleStateMatrix(rotated_rho, n, G), J_G, alpha, beta, gamma)
    appendDensityMatrixToFlatCoupledMatrix(rotated_rho, rotated_ground_rho, G, n)
    
    return rotated_rho