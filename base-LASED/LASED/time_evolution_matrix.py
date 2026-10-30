'''
This is a file to define a function to populate the time evolution matrix for a laser-atom system
Author: Manish Patel
Date created: 12/05/2021

Since v2 the matrix is assembled from a Hermitian Hamiltonian and Lindblad collapse operators by
``LASED.engine`` (vectorised, cached Wigner coefficients), which fixes the non-physical growth that
occurred for complex (elliptical) polarisation amplitudes.
'''

from LASED.symbolic_print import symbolicPrintSystem
from LASED.index import getStateLabelsFromLineNo

import numpy as np


def timeEvolutionMatrix(n, E, G, Q, Q_decay, tau, laser_wavelength, laser_intensity,
                        tau_f = None, tau_b = None, detuning = None, numeric_print = None,
                       rabi_scaling = None, rabi_factors = None, atomic_velocity = None,
                       pretty_print_eq = None, pretty_print_eq_tex = None,
                       pretty_print_eq_pdf = None, pretty_print_eq_filename = None, polarisation = None):
    """Function to create and populate the coupled differential equation matrix A for the laser-atom system.

    d(rho_flat)/dt = A rho_flat, where rho_flat[index(i, j, n)] = rho_ij.

    Returns:
        ndarray: Matrix which contains all the coefficients for the set of coupled differential equations describing a laser-atom system.
    """
    from LASED.time_evolution import buildLegacySystem  # avoid an import cycle

    system = buildLegacySystem(E, G, Q, tau, laser_wavelength, laser_intensity, tau_f, tau_b, detuning,
                               rabi_scaling, rabi_factors, polarisation=polarisation)
    A = system.liouvillian(velocity=atomic_velocity or 0.0, sparse=False)

    if numeric_print:
        for row in range(n*n):
            i, j = getStateLabelsFromLineNo(row, n)
            print("rho_dot", i, j, " = ")
            for line in np.nonzero(A[row])[0]:
                print(A[row, line], "rho", getStateLabelsFromLineNo(line, n))

    # Symbolic Printing
    if(pretty_print_eq or pretty_print_eq_pdf or pretty_print_eq_tex):
        if rabi_factors is None:
            rabi_factors = [1 for q in Q]
        symbolicPrintSystem(n, E, G, Q, Q_decay, tau_f, tau_b, detuning,
                            laser_wavelength, atomic_velocity or 0, 1 if rabi_scaling is None else rabi_scaling,
                            rabi_factors, pretty_print_eq = pretty_print_eq, pretty_print_eq_tex = pretty_print_eq_tex,
                            pretty_print_eq_pdf = pretty_print_eq_pdf, pretty_print_eq_filename = pretty_print_eq_filename)
    return A
