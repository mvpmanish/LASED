'''
This is a file to define a function to populate the time evolution matrix for a laser-atom system
Author: Manish Patel
Date created: 12/05/2021
'''

import numpy as np
from state import *
from detuning import *
from sympy import *
from sympy import Symbol
from symbolic_print import *
from half_rabi_freq import *
from index import *

def timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity,
                        tau_f = None, symbolic_print = None, numeric_print = None,
                       rabi_scaling = None, atomic_velocity = None):
    '''
    Function to create and populate the coupled differential equation matrix A for the laser-atom system.
    '''
    # Initialise matrix with zeros
    A = np.zeros((n*n,n*n), dtype = np.complex)
    
    # Calculate half-Rabi frequency
    rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength)
    if(rabi_scaling != None):
        rabi = rabi*rabi_scaling
    
    # Initialise null parameters
    if(atomic_velocity == None):
        atomic_velocity = 0
    
    # Populate A matrix
    rho_ggpp(A, n, E, G, Q, Q_decay, tau, rabi, numeric_print = numeric_print)
    rho_eepp(A, n, E, G, Q, Q_decay, tau, rabi, tau_f = tau_f, numeric_print = numeric_print)
    rho_ge(A, n, E, G, Q, Q_decay, tau, rabi, laser_wavelength, atomic_velocity, tau_f = tau_f, numeric_print = numeric_print)
    rho_eg(A, n, E, G, Q, Q_decay, tau, rabi, laser_wavelength, atomic_velocity, tau_f = tau_f, numeric_print = numeric_print)

    # Symbolic Printing
    if(symbolic_print == True):
        init_printing()
        symbolicPrintSystem(n, E, G, Q, Q_decay, tau_f, laser_wavelength, atomic_velocity)
    
    return A

def rho_ggpp(A, n, E, G, Q, Q_decay, tau, rabi, numeric_print = None):
    # rho_gg''
    for g in G:  # Start with looping over g and g'' for rho_gg''
        for gpp in G:
            row = index(g, gpp, n)  # matrix positions of rho_gg'' in 1D array
            A[row, row] += -1.j*detuning(g, gpp)  # first term in equation
            for e in E:
                column = index(g, e, n)
                for q in Q:  # Sum over all polarisations
                    A[row, column] += coupling(e, gpp, q)*1.j*rabi
            for e in E:
                column = index(e, gpp, n)
                for q in Q: 
                    A[row, column] += -1.j*coupling(e, g, q)*rabi
            for ep in E:
                for epp in E:
                    column = index(epp, ep, n)
                    column2 = index(ep, epp, n)
                    sum_decay_channels = 0
                    for gp in G:
                        for qp in Q_decay:  # Sum over decay channel polarisations
                            sum_decay_channels += coupling(epp, gp, qp)*coupling(ep, gp, qp)
                    if(sum_decay_channels != 0):
                        for qp in Q_decay:
                            A[row, column] += 1/(2*tau)*coupling(ep, gpp, qp)*coupling(epp,g, qp)/sum_decay_channels
                            A[row, column2] += 1/(2*tau)*coupling(epp, gpp, qp)*coupling(ep, g, qp)/sum_decay_channels
            if(numeric_print == True):  # Print numerical equations
                print("rho_dot", g.label, gpp.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

def rho_eepp(A, n, E, G, Q, Q_decay, tau, rabi, tau_f = None, numeric_print = None):
    # rho_ee''
    for e in E:
        for epp in E:
            row = index(e, epp, n)
            A[row, row] += -1.j*detuning(e, epp) - 1/(tau)
            if(tau_f != None):
                A[row, row] -= 1/tau_f
            for g in G:
                column = index(e, g, n)
                for q in Q:
                    A[row, column] += 1.j*coupling(epp, g, q)*rabi
            for g in G:
                column = index(g, epp, n)
                for q in Q: 
                    A[row, column] += -1.j*coupling(e, g, q)*rabi
            if(numeric_print == True):
                print("rho_dot", e.label, epp.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

def rho_ge(A, n, E, G, Q, Q_decay, tau, rabi, laser_wavelength, atomic_velocity, tau_f = None, numeric_print = None):             
    # rho_ge
    for g in G:
        for e in E: 
            row = index(g, e, n)
            A[row, row] += -1.j*dopplerDelta(e, g, w_q = angularFreq(laser_wavelength),
                                             lambda_q = laser_wavelength, v_z = atomic_velocity)  - 1/(2*tau) 
            if(tau_f != None):
                A[row, row] -= 1/(2*tau_f)
            for ep in E:
                column = index(ep, e, n)
                for q in Q: 
                    A[row, column] += -1.j*coupling(ep, g, q)*rabi
            for gp in G:
                column = index(g, gp, n)
                for q in Q: 
                    A[row, column] += 1.j*coupling(e, gp, q)*rabi
            if(numeric_print == True):
                print("rho_dot", g.label, e.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

def rho_eg(A, n, E, G, Q, Q_decay, tau, rabi, laser_wavelength, atomic_velocity, tau_f = None, numeric_print = None):              
    # rho_eg
    for e in E:
        for g in G:  
            row = index(e, g, n)
            A[row, row] += 1.j*dopplerDelta(e, g, w_q = angularFreq(laser_wavelength), 
                                            lambda_q = laser_wavelength, v_z = atomic_velocity)  - 1/(2*tau)
            if(tau_f != None):
                A[row, row] -= 1/(2*tau_f)
            for ep in E:
                column = index(e, ep, n)
                for q in Q: 
                    A[row, column] += 1.j*coupling(ep, g, q)*rabi
            for gp in G:
                column = index(gp, g, n)
                for q in Q: 
                    A[row, column] += -1.j*coupling(e, gp, q)*rabi
            if(numeric_print == True):
                print("rho_dot", e.label, g.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))