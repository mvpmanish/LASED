'''
Defines functions to symbolically print the equations of motion of the laser-atom system.
Author: Manish Patel
Date created: 07/06/2021

Since v2 the symbolic equations are generated from the same Hermitian Hamiltonian and Lindblad
collapse operators that the numerical engine uses, so the printed equations are exactly the ones solved:

    d rho/dt = -i[H, rho] + sum_q (C_q rho C_q^dag - 1/2 {C_q^dag C_q, rho}) - 1/2 {K, rho}

Symbols: Omega(e, g, q) is the (real) dipole coupling of sub-states e and g for polarisation q times the
half-Rabi frequency; the complex polarisation amplitudes appear as numbers multiplying it and are
complex-conjugated in the H_ge elements. omega_L is the laser angular frequency, delta the detuning,
k*v_z the Doppler shift, w_i the angular frequency of state i, tau the lifetime.
'''

from LASED.state import *
from LASED.detuning import *
from LASED.half_rabi_freq import coupling
from sympy import (Symbol, Function, Matrix, Eq, I, S, conjugate, sqrt, latex, expand, Float, Integer,
                   init_printing, zeros)
import numpy as np
import os

try:  # IPython is optional
    from IPython.display import display
except ImportError:  # pragma: no cover
    display = print


def _number(x):
    """A tidy sympy number from a float/complex."""
    x = complex(x)

    def part(v):
        v = round(v, 10)
        return Integer(int(v)) if v == int(v) else Float(v, 6)
    return part(x.real) + I*part(x.imag)


def symbolicEquations(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength,
                      atomic_velocity, rabi_scaling, rabi_factors):
    """Symbolic right-hand sides of the equations of motion.

    Returns:
        list of (sympy.Eq): One equation d rho_ij/dt = ... per density matrix element, ordered as index(i, j, n).
    """
    Omega = Function("Omega", real=True)
    tau = Symbol("tau", positive=True)
    omega_L, delta, k, v_z = Symbol("omega_L"), Symbol("delta"), Symbol("k"), Symbol("v_z")
    rho = Matrix(n, n, lambda i, j: Symbol(f"rho_{i+1}{j+1}"))

    H = zeros(n, n)
    for g in G:
        H[g.label-1, g.label-1] = Symbol(f"w_{g.label}")
    for e in E:
        shift = -omega_L + (k*v_z if atomic_velocity else 0) - (delta if detuning else 0)
        H[e.label-1, e.label-1] = Symbol(f"w_{e.label}") + shift
    scale = 1 if rabi_scaling is None else rabi_scaling
    factors = [1]*len(Q) if rabi_factors is None else rabi_factors
    for e in E:
        for g in G:
            V = 0
            for q, f in zip(Q, factors):
                if coupling(e, g, q) != 0:
                    V += _number(scale*f)*Omega(e.label, g.label, q)
            H[e.label-1, g.label-1] += V
            H[g.label-1, e.label-1] += conjugate(V)

    rhs = -I*(H*rho - rho*H)
    for q in Q_decay:
        Cq = zeros(n, n)
        for e in E:
            norm = sum(coupling(e, gp, qp)**2 for gp in G for qp in (-1, 0, 1))
            if norm == 0:
                continue
            for g in G:
                c = coupling(e, g, q)
                if c != 0:
                    Cq[g.label-1, e.label-1] = sqrt(_number(c**2/norm))*(1 if c > 0 else -1)/sqrt(tau)
        if any(x != 0 for x in Cq):
            CdC = Cq.H*Cq
            rhs += Cq*rho*Cq.H - (CdC*rho + rho*CdC)/2
    K = zeros(n, n)
    if tau_f is not None:
        for e in E:
            K[e.label-1, e.label-1] = 1/Symbol("tau_f", positive=True)
    if tau_b is not None:
        for g in G:
            K[g.label-1, g.label-1] = 1/Symbol("tau_b", positive=True)
    rhs += -(K*rho + rho*K)/2
    return [Eq(Symbol(f"rhodot_{i+1}{j+1}"), expand(rhs[i, j])) for i in range(n) for j in range(n)]


def symbolicPrintSystem(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength,
                        atomic_velocity, rabi_scaling, rabi_factors,
                        pretty_print_eq = None, pretty_print_eq_tex = None,
                        pretty_print_eq_pdf = None, pretty_print_eq_filename = None):
    """Prints the equations of motion of the laser-atom system in full using Sympy.

    Returns:
        list of sympy.Eq: The equations (also displayed and/or written to .tex/.pdf as requested).
    """
    if((pretty_print_eq_pdf or pretty_print_eq_tex) and (not pretty_print_eq_filename)):
        print("""Need to have a filename for the .tex/.pdf files!
        \nUse the keyword \"pretty_print_eq_filename =\" and specify your filename.
        """)
        return None
    if(pretty_print_eq):
        init_printing() # Initialise the IPython environment for printing
    equations = symbolicEquations(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength,
                                  atomic_velocity, rabi_scaling, rabi_factors)
    tex_file = None
    if(pretty_print_eq_pdf or pretty_print_eq_tex):
        tex_file = createTexFile(pretty_print_eq_filename)
    for eq in equations:
        if(pretty_print_eq):
            display(eq)
        if(tex_file):
            appendEqToTexFile(tex_file, latex(eq))
    if(tex_file):
        closeTexFile(tex_file)
        if(pretty_print_eq_pdf):
            callPdfLatex(pretty_print_eq_filename)
    return equations

def createTexFile(filename):
    """Creates a .tex file with the filename given and returns a file

    Parameters:
        filename (string): Name of the .tex file created.

    Returns:
        file: .tex file with filename
    """
    out_file = open(f'{filename}.tex', "w")
    header = "\\documentclass{article}\n\\begin{document}\n"
    header += f"\\title{{{filename}}}\n"
    header += "\\maketitle\n"
    out_file.write(header)
    return out_file

def appendEqToTexFile(file, input_str):
    """Appends an equation to a .tex file given.

    Parameters:
        file (file): file object to write to
        input_str (string): The string (equation) to be written to the .tex file
    """
    file.write("\\begin{equation}\n")
    file.write(input_str+"\n")
    file.write("\\end{equation}\n")

def closeTexFile(file):
    """Closes a .tex file which has been written to

    Parameters:
        file (file): file object to be closed
    """
    file.write("\\end{document}")
    file.close()

def callPdfLatex(filename):
    """Calls pdflatex on the system to convert the .tex file to a .pdf file. Must have pdflatex installed.

    Parameters:
        filename (string): Name of the .tex file to be converted.
    """
    os.system(f"pdflatex {filename}.tex")
