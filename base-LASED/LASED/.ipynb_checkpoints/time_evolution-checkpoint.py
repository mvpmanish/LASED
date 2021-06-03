'''
This file contains the function to calculate the evolution of the density matrix
for an atomic system interacting with a laser
'''

from constants import *
from detuning import *
from halfRabiFreq import *
from index import *
from matrix_methods import *
from state import *
import scipy.sparse.linalg as sparsela

def timeEvolution(n, E, G, Q, Q_decay, tau, tau_f, laser_power, r_sigma, n_intensity, laser_wavelength, doppler_width, doppler_detunings, time, rho0, rho_output):
    
    # Create rings of laser beam to integrate populations over
    R = np.linspace(r_sigma, 3*r_sigma, n_intensity)
    
    # Calculate doppler_spacing
    d_doppler = abs(doppler_detunings[1] - doppler_detunings[0])
    doppler_factor = d_doppler/(np.sqrt(2*PI*doppler_width*doppler_width))
    
    for doppler_delta in doppler_detunings:
        for k, r in enumerate(R, start = 0):
            # Calculate the half-rabi frequency
            rabi = halfRabiFreq(gaussianIntensity(laser_power, r_sigma, r), tau, laser_wavelength)

            # Initialise matrix with zeros
            A = np.zeros((n*n,n*n), dtype = np.complex)

            # rho_gg''
            for g in G:  # Start with looping over g and g'' for rho_gg''
                for gpp in G:
                    row = index(g, gpp, n)  # matrix positions of rho_gg'' in 1D array
                    A[row, row] += -1.j*Delta(g, gpp)  # first term in equation
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
                                for qp in Q_decay:
                                    sum_decay_channels += coupling(epp, gp, qp)*coupling(ep, gp, qp)
                            if(sum_decay_channels != 0):
                                for qp in Q_decay:
                                    A[row, column] += 1/(2*tau)*coupling(ep, gpp, qp)*coupling(epp,g, qp)/sum_decay_channels
                                    A[row, column2] += 1/(2*tau)*coupling(epp, gpp, qp)*coupling(ep, g, qp)/sum_decay_channels
        #             print("rho_dot", g.label, gpp.label, " = ")
        #             for line in range(n*n):
        #                 if (A[row, line] != 0):
        #                     print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

            # rho_ee''
            for e in E:
                for epp in E:
                    row = index(e, epp, n)
                    A[row, row] += -1.j*Delta(e, epp) - 1/(tau) - 1/(tau_f)
                    for g in G:
                        column = index(e, g, n)
                        for q in Q:
                            A[row, column] += 1.j*coupling(epp, g, q)*rabi
                            #print(A[row, column])
                    for g in G:
                        column = index(g, epp, n)
                        for q in Q: 
                            A[row, column] += -1.j*coupling(e, g, q)*rabi
                            #print(A[row, column])
        #             print("rho_dot", e.label, epp.label, " = ")
        #             for line in range(n*n):
        #                 if (A[row, line] != 0):
        #                     print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))



            # rho_ge
            for g in G:
                for e in E:
                    row = index(g, e, n)
                    A[row, row] += -1.j*doppler_delta  - 1/(2*tau) - 1/(2*tau_f)
                    for ep in E:
                        column = index(e, ep, n)
                        for q in Q: 
                            A[row, column] += -1.j*coupling(ep, g, q)*rabi
                    for gp in G:
                        column = index(g, gp, n)
                        for q in Q: 
                            A[row, column] += 1.j*coupling(e, gp, q)*rabi
        #             print("rho_dot", g.label, e.label, " = ")
        #             for line in range(n*n):
        #                 if (A[row, line] != 0):
        #                     print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

            for e in E:
                for g in G:
                    row = index(e, g, n)
                    A[row, row] += 1.j*doppler_delta - 1/(2*tau)
                    for ep in E:
                        column = index(ep, e, n)
                        for q in Q: 
                            A[row, column] += 1.j*coupling(ep, g, q)*rabi
                    for gp in G:
                        column = index(gp, g, n)
                        for q in Q: 
                            A[row, column] += -1.j*coupling(e, gp, q)*rabi
        #             print("rho_dot", e.label, g.label, " = ")
        #             for line in range(n*n):
        #                 if (A[row, line] != 0):
        #                     print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))


            # Compute the diagonalised matrix D and matrix of eigenvectors V
            D = diagonalise(A)
            V = matrixOfEigenvec(A)

            # Compute V^-1*rho(0)
            f = np.dot(la.inv(V), rho0)

            # Calculate the exponential
            for position, t in enumerate(time, start = 0):
                # Need to use expm() which computes the matrix exponential using Pade approximation
                # Differs from exp() as it useds 128 bit precision instead of 64 bit precision. This is
                # needed for a matrix exponential or loss of precision occurs
                VexpDt = np.dot(V, sparsela.expm(D*t))  # Use sparse linear algebra to compute e^Dt as this is faster
                rho_t = np.dot(VexpDt, f)

        #         if t == 0:
        #             print("At t = 0, rho(t) = ")
        #             print(rho_t)
        #         if t == (stop_time):
        #             print("At t = t_stop, rho(t) =")
        #             print(rho_t)

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
                                                                                                     