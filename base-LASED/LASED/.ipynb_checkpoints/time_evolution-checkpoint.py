'''
This file contains the function to calculate the evolution of the density matrix
for an atomic system interacting with a laser
'''

from constants import *
import numpy as np
from detuning import *
from half_rabi_freq import *
from matrix_methods import *
from time_evolution_matrix import *
import scipy.sparse.linalg as sparsela

def timeEvolution(n, E, G, Q, Q_decay, tau, laser_intensity, laser_wavelength, time, rho0, rho_output, tau_f = None, rabi_scaling = None, print_eq = None, pretty_print_eq = None, atomic_velocity = None):
    
    rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength)
            
    A = timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity, 
                        tau_f = tau_f, rabi_scaling = rabi_scaling, symbolic_print = pretty_print_eq, 
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
        
        # Append density matrix elements
        # rho(t)_ee
        for e in E:
            rho_output[index(e, e, n)][position] = abs(rho_t[index(e, e, n), 0])
        # rho(t)_gg
        for g in G:
            rho_output[index(g, g, n)][position] = abs(rho_t[index(g, g, n), 0])
        # rho(t)_eg
        for e in E:
            for g in G:
                rho_output[index(e, g, n)][position] = abs(rho_t[index(e, g, n), 0])
        #rho(t)_ge
        for g in G:
            for e in E:
                rho_output[index(g, e, n)][position] = abs(rho_t[index(e, g, n), 0])
    

def timeEvolutionGaussianAndDopplerAveraging(n, E, G, Q, Q_decay, tau, laser_power, r_sigma, n_intensity, laser_wavelength, doppler_width, doppler_detunings, time, rho0, rho_output, tau_f = None, rabi_scaling = None, print_eq = None, pretty_print_eq = None, atomic_velocity = None):
    
    if(print_eq or pretty_print_eq != None):
        print("Cannot print equations when gaussian or doppler averaging!")
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
                        tau_f = tau_f, rabi_scaling = rabi_scaling, symbolic_print = pretty_print_eq, 
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
                    rho_output[index(e, e, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(abs(rho_t[index(e, e, n), 0]))/(n_intensity*n_intensity)
                # rho(t)_gg
                for g in G:
                    rho_output[index(g, g, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(abs(rho_t[index(g, g, n), 0]))/(n_intensity*n_intensity)
                # rho(t)_eg
                for e in E:
                    for g in G:
                        rho_output[index(e, g, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(abs(rho_t[index(e, g, n), 0]))/(n_intensity*n_intensity)
                #rho(t)_ge
                for g in G:
                    for e in E:
                        rho_output[index(g, e, n)][position] += doppler_factor*np.exp(-np.power(doppler_delta/doppler_width, 2)/2)*(2*k+1)*(abs(rho_t[index(e, g, n), 0]))/(n_intensity*n_intensity)
                                                                                                     