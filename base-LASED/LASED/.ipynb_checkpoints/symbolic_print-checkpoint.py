'''
Defines functions to symbolically print the equations of motion of the laser-atom system.
Author: Manish Patel
Date created: 07/06/2021
'''

from state import *
from detuning import *
from sympy import *
from sympy import Symbol
from half_rabi_freq import *
from index import *
from decay_constant import *


def symbolicPrintSystem(n, E, G, Q, Q_decay, tau_f, detuning, laser_wavelength, atomic_velocity):
    """Prints the equations of motion of the laser-atom system in full using Sympy.
    """
    symbolicPrintRhogg(n, E, G, Q, Q_decay)
    symbolicPrintRhoee(n, E, G, Q, Q_decay, tau_f)
    symbolicPrintRhoge(n, E, G, Q, Q_decay, tau_f, detuning, laser_wavelength, atomic_velocity)
    symbolicPrintRhoeg(n, E, G, Q, Q_decay, tau_f, detuning, laser_wavelength, atomic_velocity)
    
    
def symbolicPrintRhogg(n, E, G, Q, Q_decay):
        """Prints the density matrix elements rho_gg'' for the motion of the laser-atom system using Sympy.
        """
        # rho_gg''
        for g in G:  # Start with looping over g and g'' for rho_gg''
            for gpp in G:
                rho_dot = 0
                if(delta(g, gpp) != 0):
                    rho_dot += S('Delta_{}{}*rho_{}{}'.format(g.label, gpp.label, g.label, gpp.label))
                for e in E:
                    for q in Q: 
                        if(coupling(e, gpp, q) != 0):
                            rho_dot += S('I*Omega({}, {}, {})*rho_{}{}'.format(e.label, gpp.label, q, g.label, e.label))
                for e in E:
                    column = index(e, gpp, n)
                    for q in Q: 
                        if(coupling(e, g, q) != 0):
                            rho_dot += S('-I*Omega({}, {}, {})*rho_{}{}'.format(e.label, g.label, q, e.label, gpp.label))
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
                display(Eq(S('rhodot_{}{}'.format(g.label, gpp.label)), rho_dot))

def symbolicPrintRhoee(n, E, G, Q, Q_decay, tau_f):
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
                    for q in Q:
                        if(coupling(epp, g, q) != 0):
                            rho_dot += S('I*Omega({}, {}, {})*rho_{}{}'.format(epp.label, g.label, q, e.label, epp.label))
                for g in G:
                    for q in Q: 
                        if(coupling(e, g, q) != 0):
                            rho_dot += S('-I*Omega({}, {}, {})*rho_{}{}'.format(e.label, g.label, q, g.label, epp.label))
                display(Eq(S('rhodot_{}{}'.format(e.label, epp.label)), rho_dot))
                
def symbolicPrintRhoge(n, E, G, Q, Q_decay, tau_f, detuning, laser_wavelength, atomic_velocity):
        """Prints the density matrix elements rho_ge for the motion of the laser-atom system using Sympy.
        """ 
        # rho_ge
        for g in G:
            for e in E: 
                rho_dot = S('-rho_{}{}/(2*tau)'.format(g.label, e.label))
                if(tau_f != None):
                    rho_dot += S('-rho{}{}/(2*tau_f)'.format(g.label, e.label))
                if(dopplerDelta(e, g, w_q = angularFreq(laser_wavelength), lambda_q = laser_wavelength, v_z = atomic_velocity) != 0):
                    rho_dot += S('-I*Delta({}, {}, omega_q, v_z)*rho_{}{}'.format(e.label, g.label, g.label, e.label))
                if(detuning != None):
                    rho_dot += S(f"-I*delta*rho_{g.label}{e.label}")
                for ep in E:
                    for q in Q: 
                        if(coupling(ep, g, q) != 0):
                            rho_dot += S('-I*Omega({}, {}, {})*rho_{}{}'.format(ep.label, g.label, q, ep.label, e.label))
                for gp in G:
                    for q in Q: 
                        if(coupling(e, gp, q) != 0):
                            rho_dot += S('-I*Omega({}, {}, {})*rho_{}{}'.format(e.label, gp.label, q, g.label, gp.label))
                display(Eq(S('rhodot_{}{}'.format(g.label, e.label)), rho_dot))
                
def symbolicPrintRhoeg(n, E, G, Q, Q_decay, tau_f, detuning, laser_wavelength, atomic_velocity):
        """Prints the density matrix elements rho_eg for the motion of the laser-atom system using Sympy.
        """ 
        # rho_eg
        for e in E:
            for g in G:
                rho_dot = S('-rho_{}{}/(2*tau)'.format(e.label, g.label))
                if(tau_f != None):
                    rho_dot += S('-rho{}{}/(2*tau_f)'.format(e.label, g.label))
                if(dopplerDelta(e, g, w_q = angularFreq(laser_wavelength), lambda_q = laser_wavelength, v_z = atomic_velocity) != 0):
                    rho_dot += S('I*Delta({}, {}, omega_q, v_z)*rho_{}{}'.format(e.label, g.label, e.label, g.label))
                if(detuning != None):
                    rho_dot += S(f"I*delta*rho_{e.label}{g.label}")
                for ep in E:
                    for q in Q: 
                        if(coupling(ep, g, q) != 0):
                            rho_dot += S('I*Omega({}, {}, {})*rho_{}{}'.format(ep.label, g.label, q, e.label, ep.label))
                for gp in G:
                    for q in Q: 
                        if(coupling(e, g, q) != 0):
                            rho_dot += S('I*Omega({}, {}, {})*rho_{}{}'.format(e.label, g.label, q, gp.label, g.label))
                display(Eq(S('rhodot_{}{}'.format(e.label, g.label)), rho_dot))