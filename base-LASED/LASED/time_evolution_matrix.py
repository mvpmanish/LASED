'''
This is a file to define a function to populate the time evolution matrix for a laser-atom system
Author: Manish Patel
Date created: 12/05/2021
'''

from LASED.state import *
from LASED.detuning import *
from LASED.symbolic_print import *
from LASED.half_rabi_freq import *
from LASED.decay_constant import *
from LASED.index import *

from sympy import *
from sympy import Symbol
import numpy as np

def timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity,
                        tau_f = None, detuning = None, symbolic_print = None, numeric_print = None,
                       rabi_scaling = None, rabi_factors = None, atomic_velocity = None):
    """Function to create and populate the coupled differential equation matrix A for the laser-atom system.
    
    Returns:
        (ndarray) Matrix which contains all thera coefficients for the set of coupled differential equations describing a laser-atom system. 
    """
    
    # Initialise matrix with zeros
    A = np.zeros((n*n,n*n), dtype = np.complex)
    
    # Calculate half-Rabi frequency
    rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength)
    if(rabi_scaling != None):  # For normalising the rabi frequency
        rabi = rabi*rabi_scaling
    else:
        rabi_scaling = 1  # For symbolic printing
    
    if(rabi_factors != None):
        if(len(rabi_factors) != len(Q)):
            print("rabi_factors must be the same length as Q! Each element of Q is multiplied by the corresponding rabi_factor.")
    else:
        rabi_factors = [1 for q in Q]  # Set Rabi factors to 1
    
    # Initialise null parameters
    if(atomic_velocity == None):
        atomic_velocity = 0
    
    # Populate A matrix
    rho_ggpp(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, numeric_print = numeric_print)
    rho_eepp(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, tau_f = tau_f, numeric_print = numeric_print)
    rho_ge(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, laser_wavelength, atomic_velocity, tau_f = tau_f, detuning = detuning, numeric_print = numeric_print)
    rho_eg(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, laser_wavelength, atomic_velocity, tau_f = tau_f, detuning = detuning, numeric_print = numeric_print)

    # Symbolic Printing
    if(symbolic_print == True):
        init_printing()
        symbolicPrintSystem(n, E, G, Q, Q_decay, tau_f, detuning, 
                            laser_wavelength, atomic_velocity, rabi_scaling, rabi_factors)
    
    return A

def rho_ggpp(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, numeric_print = None):
    """ Function to populate the matrix A with coefficients for populations and atomic coherences of the ground states. 
    """
    # rho_gg''
    for g in G:  # Start with looping over g and g'' for rho_gg''
        for gpp in G:
            row = index(g, gpp, n)  # matrix positions of rho_gg'' in 1D array
            A[row, row] += -1.j*delta(g, gpp)  # first term in equation
            for e in E:
                column = index(g, e, n)
                for i,q in enumerate(Q):  # Sum over all polarisations
                    A[row, column] += coupling(e, gpp, q)*1.j*rabi*rabi_factors[i]
            for e in E:
                column = index(e, gpp, n)
                for i,q in enumerate(Q): 
                    A[row, column] += -1.j*coupling(e, g, q)*rabi*rabi_factors[i]
            for ep in E:
                for epp in E:
                    column = index(epp, ep, n)
                    column2 = index(ep, epp, n)
                    sum_decay_channels = 0
                    for gp in G:
                        for qp in Q_decay:  # Sum over decay channel polarisations
                            sum_decay_channels += abs(coupling(epp, gp, qp)*coupling(ep, gp, qp))
                    if(sum_decay_channels != 0):
                        if(ep.label == epp.label):
                            for qp in Q_decay:
                                A[row, column] += 1/(2*tau)*abs(coupling(ep, gpp, qp)*coupling(epp,g, qp))/sum_decay_channels
                                A[row, column2] += 1/(2*tau)*abs(coupling(epp, gpp, qp)*coupling(ep, g, qp))/sum_decay_channels
                        else:  
                        # Then this is a vertical coherence and the generalised decay constant must be evaluated
                            decay_const = generalisedDecayConstant(ep, epp, gpp, G, Q_decay)/(2*tau)  # Divide by two to take into account double counting of e'e'' and e''e'
                            A[row, column] += decay_const
                            A[row, column2] += decay_const
                        
            if(numeric_print == True):  # Print numerical equations
                print("rho_dot", g.label, gpp.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

def rho_eepp(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, tau_f = None, numeric_print = None):
    """ Function to populate the matrix A with coefficients for populations and atomic coherences of the excited states. 
    """
    # rho_ee''
    for e in E:
        for epp in E:
            row = index(e, epp, n)
            A[row, row] += -1.j*delta(e, epp) - 1/(tau)
            if(tau_f != None):
                A[row, row] -= 1/tau_f
            for g in G:
                column = index(e, g, n)
                for i,q in enumerate(Q):
                    A[row, column] += 1.j*coupling(epp, g, q)*rabi*rabi_factors[i]
            for g in G:
                column = index(g, epp, n)
                for i,q in enumerate(Q): 
                    A[row, column] += -1.j*coupling(e, g, q)*rabi*rabi_factors[i]
            if(numeric_print == True):
                print("rho_dot", e.label, epp.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

def rho_ge(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, laser_wavelength, atomic_velocity, tau_f = None, detuning = None, numeric_print = None):
    """ Function to populate the matrix A with coefficients for optical coherences between ground and excited states. 
    """
    # rho_ge
    for g in G:
        for e in E: 
            row = index(g, e, n)
            A[row, row] += -1.j*dopplerDelta(e, g, w_q = angularFreq(laser_wavelength),
                                             lambda_q = laser_wavelength, v_z = atomic_velocity)  - 1/(2*tau)
            if(detuning != None):
                A[row, row] += -1.j*detuning
            if(tau_f != None):
                A[row, row] -= 1/(2*tau_f)
            for ep in E:
                column = index(ep, e, n)
                for i,q in enumerate(Q): 
                    A[row, column] += -1.j*coupling(ep, g, q)*rabi*rabi_factors[i]
            for gp in G:
                column = index(g, gp, n)
                for i,q in enumerate(Q): 
                    A[row, column] += 1.j*coupling(e, gp, q)*rabi*rabi_factors[i]
            if(numeric_print == True):
                print("rho_dot", g.label, e.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

def rho_eg(A, n, E, G, Q, Q_decay, tau, rabi, rabi_factors, laser_wavelength, atomic_velocity, tau_f = None, detuning = None, numeric_print = None):
    """ Function to populate the matrix A with coefficients for optical coherences between excited and ground states. 
    """
    # rho_eg
    for e in E:
        for g in G:  
            row = index(e, g, n)
            A[row, row] += 1.j*dopplerDelta(e, g, w_q = angularFreq(laser_wavelength), 
                                            lambda_q = laser_wavelength, v_z = atomic_velocity)  - 1/(2*tau)
            if(detuning != None):
                A[row, row] += 1.j*detuning
            if(tau_f != None):
                A[row, row] -= 1/(2*tau_f)
            for ep in E:
                column = index(e, ep, n)
                for i,q in enumerate(Q): 
                    A[row, column] += 1.j*coupling(ep, g, q)*rabi*rabi_factors[i]
            for gp in G:
                column = index(gp, g, n)
                for i,q in enumerate(Q): 
                    A[row, column] += -1.j*coupling(e, gp, q)*rabi*rabi_factors[i]
            if(numeric_print == True):
                print("rho_dot", e.label, g.label, " = ")
                for line in range(n*n):
                    if (A[row, line] != 0):
                        print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))