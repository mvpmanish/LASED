b'''
Defines functions to symbolically print the equations of motion of the laser-atom system.
Author: Manish Patel
Date created: 07/06/2021
'''

from LASED.state import *
from LASED.detuning import *
from sympy import *
from sympy import Symbol
from LASED.half_rabi_freq import *
from LASED.index import *
from LASED.decay_constant import *
import os


def symbolicPrintSystem(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength,
                        atomic_velocity, rabi_scaling, rabi_factors,
                        pretty_print_eq = None, pretty_print_eq_tex = None,
                        pretty_print_eq_pdf = None, pretty_print_eq_filename = None):
    """Prints the equations of motion of the laser-atom system in full using Sympy.
    """
    if(pretty_print_eq):
        init_printing() # Initialise the IPython environment for printing
    if((pretty_print_eq_pdf or pretty_print_eq_tex) and (not pretty_print_eq_filename)):
        # Catch the case where the filename is not give
        print("""Need to have a filename for the .tex/.pdf files!
        \nUse the keyword \"pretty_print_eq_filename =\" and specify your filename.
        """)
    else:
        tex_file = None  # Placeholder for a file object
        tex_file_handling = pretty_print_eq_pdf or pretty_print_eq_tex # bool to handle file actions
        if(tex_file_handling):
            # Create the tex file with the user's defined file name
            tex_file = createTexFile(pretty_print_eq_filename)
        symbolicPrintRhogg(n, E, G, Q, Q_decay, tau_b, rabi_scaling, rabi_factors,
                            pretty_print_eq = pretty_print_eq, pretty_print_eq_file = tex_file)
        symbolicPrintRhoee(n, E, G, Q, Q_decay, tau_f, rabi_scaling, rabi_factors,
                            pretty_print_eq = pretty_print_eq, pretty_print_eq_file = tex_file)
        symbolicPrintRhoge(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength,
                            atomic_velocity, rabi_scaling, rabi_factors,
                            pretty_print_eq = pretty_print_eq, pretty_print_eq_file = tex_file)
        symbolicPrintRhoeg(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength,
                            atomic_velocity, rabi_scaling, rabi_factors,
                            pretty_print_eq = pretty_print_eq, pretty_print_eq_file = tex_file)
        if(tex_file_handling):
            closeTexFile(tex_file)
            if(pretty_print_eq_pdf):
                callPdfLatex(pretty_print_eq_filename)

def symbolicPrintRhogg(n, E, G, Q, Q_decay, tau_b, rabi_scaling, rabi_factors,
                        pretty_print_eq = None, pretty_print_eq_file = None):
        """Prints the density matrix elements rho_gg'' for the motion of the laser-atom system using Sympy.
        """
        # rho_gg''
        for g in G:  # Start with looping over g and g'' for rho_gg''
            for gpp in G:
                rho_dot = 0
                if(tau_b != None):
                    rho_dot += S('-rho{}{}/tau_b'.format(g.label, gpp.label))
                if(delta(g, gpp) != 0):
                    rho_dot += S('Delta_{}{}*rho_{}{}'.format(g.label, gpp.label, g.label, gpp.label))
                for e in E:
                    for i,q in enumerate(Q):
                        if(coupling(e, gpp, q) != 0):
                            rho_dot += S('I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i], e.label, gpp.label, q, g.label, e.label))
                for e in E:
                    column = index(e, gpp, n)
                    for i,q in enumerate(Q):
                        if(coupling(e, g, q) != 0):
                            rho_dot += S('-I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i],e.label, g.label, q, e.label, gpp.label))
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
                                    rho_dot += S('{}/(2*tau)*rho_{}{} + {}/(2*tau)*rho_{}{}'.format(abs(coupling(ep, gpp, qp)*coupling(epp,g, qp))/sum_decay_channels,epp.label, ep.label,
                                                                                 abs(coupling(epp, gpp, qp)*coupling(ep, g, qp))/sum_decay_channels, ep.label, epp.label))
                            else:
                                # Gerneralised decay constant
                                rho_dot += S(f'({generalisedDecayConstant(ep, epp, gpp, G, Q_decay)}/(2*tau))*rho_{epp.label}{ep.label} + ({generalisedDecayConstant(ep, epp, gpp, G, Q_decay)}/(2*tau))*rho_{ep.label}{epp.label}')
                eq = Eq(S('rhodot_{}{}'.format(g.label, gpp.label)), rho_dot)
                if(pretty_print_eq):
                    display(eq)
                if(pretty_print_eq_file):
                     appendEqToTexFile(pretty_print_eq_file, latex(eq))


def symbolicPrintRhoee(n, E, G, Q, Q_decay, tau_f, rabi_scaling, rabi_factors,
                        pretty_print_eq = None, pretty_print_eq_file = None):
        """Prints the density matrix elements rho_ee'' for the motion of the laser-atom system using Sympy.
        """
        # rho_ee''
        for e in E:
            for epp in E:
                rho_dot = S('-1/tau*rho_{}{}'.format(e.label, epp.label))
                if(tau_f != None):
                    rho_dot += S('-rho{}{}/tau_f'.format(e.label, epp.label))
                if(delta(e, epp) != 0):
                    rho_dot += S('Delta_{}{}*rho_{}{}'.format(e.label, epp.label, e.label, epp.label))
                for g in G:
                    for i,q in enumerate(Q):
                        if(coupling(epp, g, q) != 0):
                            rho_dot += S('I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i], epp.label, g.label, q, e.label, g.label))
                for g in G:
                    for i,q in enumerate(Q):
                        if(coupling(e, g, q) != 0):
                            rho_dot += S('-I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i], e.label, g.label, q, g.label, epp.label))
                eq = Eq(S('rhodot_{}{}'.format(e.label, epp.label)), rho_dot)
                if(pretty_print_eq):
                    display(eq)
                if(pretty_print_eq_file):
                     appendEqToTexFile(pretty_print_eq_file, latex(eq))

def symbolicPrintRhoge(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength, atomic_velocity,
                        rabi_scaling, rabi_factors, pretty_print_eq = None, pretty_print_eq_file = None):
        """Prints the density matrix elements rho_ge for the motion of the laser-atom system using Sympy.
        """
        # rho_ge
        for g in G:
            for e in E:
                rho_dot = S('-rho_{}{}/(2*tau)'.format(g.label, e.label))
                if(tau_f != None):
                    rho_dot += S('-rho{}{}/(2*tau_f)'.format(g.label, e.label))
                if(tau_b != None):
                    rho_dot += S('-rho{}{}/(2*tau_b)'.format(g.label, e.label))
                if(dopplerDelta(e, g, w_q = angularFreq(laser_wavelength), lambda_q = laser_wavelength, v_z = atomic_velocity) != 0):
                    rho_dot += S('-I*Delta({}, {}, omega_q, v_z)*rho_{}{}'.format(e.label, g.label, g.label, e.label))
                if(detuning != None):
                    rho_dot += S(f"-I*delta*rho_{g.label}{e.label}")
                for ep in E:
                    for i,q in enumerate(Q):
                        if(coupling(ep, g, q) != 0):
                            rho_dot += S('-I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i], ep.label, g.label, q, ep.label, e.label))
                for gp in G:
                    for i,q in enumerate(Q):
                        if(coupling(e, gp, q) != 0):
                            rho_dot += S('-I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i], e.label, gp.label, q, g.label, gp.label))
                eq = Eq(S('rhodot_{}{}'.format(g.label, e.label)), rho_dot)
                if(pretty_print_eq):
                    display(eq)
                if(pretty_print_eq_file):
                     appendEqToTexFile(pretty_print_eq_file, latex(eq))

def symbolicPrintRhoeg(n, E, G, Q, Q_decay, tau_f, tau_b, detuning, laser_wavelength, atomic_velocity,
                        rabi_scaling, rabi_factors, pretty_print_eq = None, pretty_print_eq_file = None):
        """Prints the density matrix elements rho_eg for the motion of the laser-atom system using Sympy.
        """
        # rho_eg
        for e in E:
            for g in G:
                rho_dot = S('-rho_{}{}/(2*tau)'.format(e.label, g.label))
                if(tau_f != None):
                    rho_dot += S('-rho{}{}/(2*tau_f)'.format(e.label, g.label))
                if(tau_b != None):
                    rho_dot += S('-rho{}{}/(2*tau_b)'.format(e.label, g.label))
                if(dopplerDelta(e, g, w_q = angularFreq(laser_wavelength), lambda_q = laser_wavelength, v_z = atomic_velocity) != 0):
                    rho_dot += S('I*Delta({}, {}, omega_q, v_z)*rho_{}{}'.format(e.label, g.label, e.label, g.label))
                if(detuning != None):
                    rho_dot += S(f"I*delta*rho_{e.label}{g.label}")
                for ep in E:
                    for i,q in enumerate(Q):
                        if(coupling(ep, g, q) != 0):
                            rho_dot += S('I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i], ep.label, g.label, q, e.label, ep.label))
                for gp in G:
                    for q in Q:
                        if(coupling(e, g, q) != 0):
                            rho_dot += S('I*{}*{}*Omega({}, {}, {})*rho_{}{}'.format(rabi_scaling, rabi_factors[i], e.label, g.label, q, gp.label, g.label))
                eq = Eq(S('rhodot_{}{}'.format(e.label, g.label)), rho_dot)
                if(pretty_print_eq):
                    display(eq)
                if(pretty_print_eq_file):
                     appendEqToTexFile(pretty_print_eq_file, latex(eq))

def createTexFile(filename):
    """Creates a .tex file with the filename given and returns a file

    Parameters:
        filename (string): Name of the .tex file created.

    Returns:
        file: .tex file with filename
    """
    out_file = open(f'{filename}.tex', "w")
    header = "\documentclass{article}\n\\begin{document}\n"
    header += f"\\title{{{filename}}}\n"
    header += "\maketitle\n"
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
    file.write("\end{equation}\n")

def closeTexFile(file):
    """Closes a .tex file which has been written to

    Parameters:
        file (file): file object to be closed
    """
    file.write("\end{document}")
    file.close()

def callPdfLatex(filename):
    """Calls pdflatex on the system to convert the .tex file to a .pdf file. Must have pdflatex installed.

    Parameters:
        filename (string): Name of the .tex file to be converted.
    """
    os.system(f"pdflatex {filename}.tex")
