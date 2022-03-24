'''Define functions for rotating density matrices in the QED simulation of a laser-atom system.
Author: Manish Patel
Date created: 12/05/2021
'''

from LASED.density_matrix import *
import numpy as np
import math

def wigner_D(J, mu, m, alpha, beta, gamma):
    """Calculates the Wigner D-matrix for rotation by Eueler angles (alpha, beta, gamma).
    
    Parameters:
        J (int): total angular momentum quantum number of the state which will be rotated with the 
        resulting D-matrix.
        alpha (float): rotation around z-axis in radians.
        beta (float): rotation about the y'-axis in radians.
        gamma (float): rotation about the z''-axis in radians.
    
    Returns:
        ndarray: A square matrix of size 2J+1.
    """
    alpha_const = math.cos(-mu*alpha)+1.j*math.sin(-mu*alpha)
    gamma_const = math.cos(-m*gamma)+1.j*math.sin(-m*gamma)
    D = alpha_const*small_Wigner_D(J, beta, mu, m)*gamma_const
    
    return D

def small_Wigner_D(J, beta, mp, m):
    """Calculates the small Wigner D-matrix elements for rotation.
    
    Parameters:
        J (int): total angular momentum quantum number of the state which will be rotated with the resulting D-matrix.
        beta (float): rotation about the y'-axis in radians.
        mp (int): row number of element in the Wigner D-matrix
        m (int): column number oif element in the Wigner D-matrix
    
    Returns:
        float: The small Wigner D-matrix elements
    """
    const = np.sqrt((math.factorial(J+mp))*math.factorial(J-mp)*math.factorial(J+m)*math.factorial(J-m))
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

def createDictionaryOfSubStates(E, G):
    """Creates a dictionary of sub-states with (F, m) as the key and the State object as the value.
    
    Parameters:
        E (list of States): List of excited sub-states in the laser-atom system.
        G (list of States): List of ground sub-states in the laser-atom system
    
    Returns:
        dictionary: A dictionary with (F, m) as the keys and the corresponding State object as the value
    """
    sub_state_list = []  # List of values
    key_list = []
    for e in E:  # Appen excited states
        sub_state_list.append(e)
    for g in G:  # Append ground states
        sub_state_list.append(g)
    for sub_state in sub_state_list:  # Append F and m values to create key values
        key_list.append((sub_state.F, sub_state.m))  # Each key is a tuple of (F, m)
    return dict(zip(key_list, sub_state_list))

def rotateElement(rho, i, j, n, sub_state_dict, alpha, beta, gamma):
    """Rotates an element of a density matrix rho_ij by the euler angles given.
    
    Parameters:
        rho (complex): The density matrix to be rotated.
        i (state): A sub-state of the laser-atom system.
        j (state): A sub-state of the laser-atom system.
        sub_state_dict (dict of State): Dictionary of States with (F,m) tuple as the keys and the corresponding State object as the value.
        n (int): Total number of sub-states in the system.
        alpha (float): rotation around z-axis in radians.
        beta (float): rotation about the y'-axis in radians.
        gamma (float): rotation about the z''-axis in radians.
    
    Returns:
        complex: The rotated density matrix element.
    """
    F = i.F
    Fp = j.F
    m = i.m
    mp = j.m
    mu_list = np.linspace(-F, F, 2*F+1, dtype = int)
    mup_list = np.linspace(-Fp, Fp, 2*Fp+1, dtype = int)
    rho_new_frame = 0  # Initialise the summation to zero
    for mu in mu_list:
        for mup in mup_list:
            D_mu_m = wigner_D(F, mu, m, alpha, beta, gamma)
            D_mup_mp = wigner_D(Fp, mup, mp, alpha, beta, gamma)
            rho_FmuFmup = rho[index(sub_state_dict.get((F, mu)), sub_state_dict.get((Fp, mup)), n), 0]  # Retrieve the density matrix element
            rho_new_frame += np.conj(D_mu_m)*rho_FmuFmup*D_mup_mp          
    return rho_new_frame

def rotateFlatDensityMatrix(flat_rho, n, E, G, alpha, beta, gamma):
    """Rotate the excited and ground state populations by the Euler angles.
    
    Parameters:
        flat_rho (list): A flattened 2D density matrix
        n (int): Number of substates which compose the density matrix flat_rho
        E (list of States): list of excited State objects which compose flat_rho
        G (list of States): list of ground State objects which compose flat_rho
        alpha (float): rotation around z-axis in radians.
        beta (float): rotation about the y'-axis in radians.
        gamma (float): rotation about the z''-axis in radians.
    
    Returns:
        list of lists: A rotated flattened 2D density matrix
    """
    rotated_rho = np.zeros((n*n, 1), dtype = complex)  # Placeholder for rotated density matrix
    sub_state_dict = createDictionaryOfSubStates(E, G)  # Dictionary to retrive states by only using (F, m) values
    # Rotate the excited state populations and atomic coherences
    # rho_gg''
    for g in G:
        for gpp in G:
            rho_ggpp_new_frame = rotateElement(flat_rho, g, gpp, n, sub_state_dict, alpha, beta, gamma)
            rotated_rho[index(g, gpp, n), 0] = rho_ggpp_new_frame
    # rho_ee''
    for e in E:
        for epp in E:
            rho_eepp_new_frame = rotateElement(flat_rho, e, epp, n, sub_state_dict, alpha, beta, gamma)
            rotated_rho[index(e, epp, n), 0] = rho_eepp_new_frame
    # rho_ge
    for g in G:
        for e in E:
            rho_ge_new_frame = rotateElement(flat_rho, g, e, n, sub_state_dict, alpha, beta, gamma)
            rotated_rho[index(g, e, n), 0] = rho_ge_new_frame
    # rho_eg
    for e in E:
        for g in G:
            rho_eg_new_frame = rotateElement(flat_rho, e, g, n, sub_state_dict, alpha, beta, gamma)
            rotated_rho[index(e, g, n), 0] = rho_eg_new_frame

    return rotated_rho