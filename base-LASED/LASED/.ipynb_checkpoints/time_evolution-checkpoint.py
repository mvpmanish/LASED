'''
This file contains the function to calculate the time evolution of the density matrix
for an atomic system interacting with a laser.
'''

from LASED.constants import *
from LASED.detuning import *
from LASED.half_rabi_freq import *
from LASED.matrix_methods import *
from LASED.time_evolution_matrix import *

import numpy as np
import scipy.linalg as la

def timeEvolution(n, E, G, Q, Q_decay, tau, laser_intensity, laser_wavelength, time, rho0, rho_output, tau_f = None, detuning = None, rabi_scaling = None, rabi_factors = None, print_eq = None, pretty_print_eq = None, atomic_velocity = None):
    """Calculates the time evolution of a laser-atom system.
    
    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output.
    """
    
    rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength)
            
    A = timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity, 
                        tau_f = tau_f, detuning = detuning, rabi_scaling = rabi_scaling,
                        rabi_factors = rabi_factors, symbolic_print = pretty_print_eq, 
                        numeric_print = print_eq, atomic_velocity = atomic_velocity)

    # Compute the diagonalised matrix D and matrix of eigenvectors V
    D = diagonalise(A)
    V = matrixOfEigenvec(A)
    f = np.dot(la.inv(V), rho0)  # Compute V^-1*rho(0)
    
    # Calculate the exponential
    for position, t in enumerate(time, start = 0):
        # Use expm() which computes the matrix exponential using the Pade approximation
        expS = la.expm(D*t) # Compute exp(D*t), this takes the longest so needs speeding up
        VexpDt = np.dot(V, expS)  
        rho_t = np.dot(VexpDt, f)
        
        # Append density matrix elements
        # rho(t)_ee
        for e in E:
            for ep in E:
                rho_output[index(e, ep, n)][position] = rho_t[index(e, ep, n), 0]
        # rho(t)_gg
        for g in G:
            for gp in G:
                rho_output[index(g, gp, n)][position] = rho_t[index(g, gp, n), 0]
        # rho(t)_eg
        for e in E:
            for g in G:
                rho_output[index(e, g, n)][position] = rho_t[index(e, g, n), 0]
        #rho(t)_ge
        for g in G:
            for e in E:
                rho_output[index(g, e, n)][position] = rho_t[index(g, e, n), 0]

def timeEvolutionDopplerAveraging(n, E, G, Q, Q_decay, tau, laser_intensity, laser_wavelength, doppler_width, doppler_detunings, time, rho0, rho_output, tau_f = None, detuning = None, rabi_scaling = None, rabi_factors = None, print_eq = None, pretty_print_eq = None, atomic_velocity = None):
    """Calculates the time evolution of a laser-atom system with a Gaussian doppler profile for the atoms.
    
    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output. 
    """
    
    if(print_eq or pretty_print_eq != None):
        print("Cannot print equations when beam profile or doppler averaging!")
        print_eq = None
        pretty_print_eq = None
    
    # Calculate doppler_spacing
    d_doppler = abs(doppler_detunings[1] - doppler_detunings[0])
    doppler_factor = d_doppler/(np.sqrt(2*PI*doppler_width*doppler_width))
    
    for doppler_delta in doppler_detunings:
            # Calculate the half-rabi frequency
            rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength)
            
            A = timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity, 
                        tau_f = tau_f, detuning = detuning, rabi_scaling = rabi_scaling, 
                        rabi_factors = rabi_factors, symbolic_print = pretty_print_eq, 
                        numeric_print = print_eq, atomic_velocity = atomic_velocity)

            # Compute the diagonalised matrix D and matrix of eigenvectors V
            D = diagonalise(A)
            V = matrixOfEigenvec(A)
            f = np.dot(la.inv(V), rho0)  # Compute V^-1*rho(0)

            # Calculate the exponential
            for position, t in enumerate(time, start = 0):
                # Use expm() which computes the matrix exponential using the Pade approximation
                VexpDt = np.dot(V, sparsela.expm(D*t))  # Use sparse linear algebra to compute e^Dt as this is faster
                rho_t = np.dot(VexpDt, f)

                # Append density matrix for each ring and each detuning fraction
                # rho(t)_ee
                for e in E:
                    for ep in E:
                        rho_output[index(e, ep, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(rho_t[index(e, ep, n), 0])
                # rho(t)_gg
                for g in G:
                    for gp in G:
                        rho_output[index(g, gp, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(rho_t[index(g, gp, n), 0])
                # rho(t)_eg
                for e in E:
                    for g in G:
                        rho_output[index(e, g, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(rho_t[index(e, g, n), 0])
                #rho(t)_ge
                for g in G:
                    for e in E:
                        rho_output[index(g, e, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(rho_t[index(e, g, n), 0])
                        
def timeEvolutionGaussianAveraging(n, E, G, Q, Q_decay, tau, laser_power, r_sigma, n_intensity, laser_wavelength, time, rho0, rho_output, tau_f = None, detuning = None, rabi_scaling = None, rabi_factors = None, print_eq = None, pretty_print_eq = None, atomic_velocity = None):
    """Calculates the time evolution of a laser-atom system with a Gaussian laser beam profile.
    
    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output. 
    """
    
    if(print_eq or pretty_print_eq != None):
        print("Cannot print equations when beam profile or doppler averaging!")
        print_eq = None
        pretty_print_eq = None
    
    # Create rings of laser beam to integrate populations over
    R = np.linspace(r_sigma, 3*r_sigma, n_intensity)
    
    for k, r in enumerate(R, start = 0):
            # Calculate the half-rabi frequency
            laser_intensity = gaussianIntensity(laser_power, r_sigma, r)
            rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength)
            
            A = timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity, 
                        tau_f = tau_f, rabi_scaling = rabi_scaling, rabi_factors = rabi_factors, 
                        detuning = detuning, symbolic_print = pretty_print_eq, 
                        numeric_print = print_eq, atomic_velocity = atomic_velocity)

            # Compute the diagonalised matrix D and matrix of eigenvectors V
            D = diagonalise(A)
            V = matrixOfEigenvec(A)
            f = np.dot(la.inv(V), rho0)  # Compute V^-1*rho(0)

            # Calculate the exponential
            for position, t in enumerate(time, start = 0):
                # Use expm() which computes the matrix exponential using the Pade approximation
                VexpDt = np.dot(V, sparsela.expm(D*t))  # Use sparse linear algebra to compute e^Dt as this is faster
                rho_t = np.dot(VexpDt, f)

                # Append density matrix for each ring
                # rho(t)_ee
                for e in E:
                    for ep in E:
                        rho_output[index(e, ep, n)][position] += (2*k+1)*(rho_t[index(e, ep, n), 0])/(n_intensity*n_intensity)
                # rho(t)_gg
                for g in G:
                    for gp in G:
                        rho_output[index(g, gp, n)][position] += (2*k+1)*(rho_t[index(g, gp, n), 0])/(n_intensity*n_intensity)
                # rho(t)_eg
                for e in E:
                    for g in G:
                        rho_output[index(e, g, n)][position] += (2*k+1)*(rho_t[index(e, g, n), 0])/(n_intensity*n_intensity)
                #rho(t)_ge
                for g in G:
                    for e in E:
                        rho_output[index(g, e, n)][position] += (2*k+1)*(rho_t[index(e, g, n), 0])/(n_intensity*n_intensity)
    

def timeEvolutionGaussianAndDopplerAveraging(n, E, G, Q, Q_decay, tau, laser_power, r_sigma, n_intensity, laser_wavelength, doppler_width, doppler_detunings, time, rho0, rho_output, tau_f = None, detuning = None, rabi_scaling = None, rabi_factors = None, print_eq = None, pretty_print_eq = None, atomic_velocity = None):
    """Calculates the time evolution of a laser-atom system with a Gaussian doppler profile for the atoms and a Gaussian laser beam profile.
    
    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output. 
    """
    
    if(print_eq or pretty_print_eq != None):
        print("Cannot print equations when beam profile or doppler averaging!")
        print_eq = None
        pretty_print_eq = None
    
    # Create rings of laser beam to integrate populations over
    R = np.linspace(r_sigma, 3*r_sigma, n_intensity)
    
    # Calculate doppler_spacing
    d_doppler = abs(doppler_detunings[1] - doppler_detunings[0])
    doppler_factor = d_doppler/(np.sqrt(2*PI*doppler_width*doppler_width))
    
    for doppler_delta in doppler_detunings:
        for k, r in enumerate(R, start = 0):
            # Calculate the half-rabi frequency
            laser_intensity = gaussianIntensity(laser_power, r_sigma, r)
            rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength)
            
            A = timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity, 
                        tau_f = tau_f, rabi_scaling = rabi_scaling, rabi_factors = rabi_factors, 
                        detuning = detuning, symbolic_print = pretty_print_eq, 
                        numeric_print = print_eq, atomic_velocity = atomic_velocity)

            # Compute the diagonalised matrix D and matrix of eigenvectors V
            D = diagonalise(A)
            V = matrixOfEigenvec(A)
            f = np.dot(la.inv(V), rho0)  # Compute V^-1*rho(0)

            # Calculate the exponential
            for position, t in enumerate(time, start = 0):
                # Use expm() which computes the matrix exponential using the Pade approximation
                VexpDt = np.dot(V, sparsela.expm(D*t))  # Use sparse linear algebra to compute e^Dt as this is faster
                rho_t = np.dot(VexpDt, f)

                # Append density matrix for each ring and each detuning fraction
                # rho(t)_ee
                for e in E:
                    for ep in E:
                        rho_output[index(e, ep, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(rho_t[index(e, ep, n), 0])/(n_intensity*n_intensity)
                # rho(t)_gg
                for g in G:
                    for gp in G:
                        rho_output[index(g, gp, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(rho_t[index(g, gp, n), 0])/(n_intensity*n_intensity)
                # rho(t)_eg
                for e in E:
                    for g in G:
                        rho_output[index(e, g, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(rho_t[index(e, g, n), 0])/(n_intensity*n_intensity)
                #rho(t)_ge
                for g in G:
                    for e in E:
                        rho_output[index(g, e, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(rho_t[index(e, g, n), 0])/(n_intensity*n_intensity)
                                                                                                     