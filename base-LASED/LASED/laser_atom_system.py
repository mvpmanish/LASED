'''
The LaserAtomSystem class definition.
'''

import time_evolution as te
import numpy as np
import index as ix

class LaserAtomSystem:
    """A user-defined laser field acting on a user-defined atomic system
    """
    
    # Class variables
    Q_decay = [1, 0 ,-1]
    rho_t = []
    
    def __init__(self, E, G, tau, Q, laser_wavelength, laser_intensity = None, 
                 laser_power = None, tau_f = None):
        self.E = E
        self.G = G
        self.tau = tau
        self.Q = Q
        self.laser_wavelength = laser_wavelength
        self.tau_f = tau_f
        self.laser_intensity = laser_intensity
        self.laser_power = laser_power
        self.rho0 = np.zeros((self.n*self.n, 1), dtype = complex)
        
    @property
    def n(self):
        return int(len(self.G)+len(self.E))
    
    
    
    def setRho0(self, i, j, value):
        if(value > 1):
            print("Cannot set an element of a density matrix > 1!")
            return
        else:
            row = ix.index(i, j, self.n)
            self.rho0[row, 0] = value
    
    def timeEvolution(self, time, beam_profile_averaging = None, doppler_averaging = None, 
                     pretty_print_eq = None, print_eq = None, rabi_scaling = None,
                     atomic_velocity = None, r_sigma = None, n_beam_averaging = None, 
                     doppler_width = None, doppler_detunings = None):
        """ Evolves the laser-atom system over time.
        """
        n = self.n
        E = self.E
        G = self.G
        Q = self.Q
        tau = self.tau
        laser_power = self.laser_power
        laser_intensity = self.laser_intensity
        laser_wavelength = self.laser_wavelength
        rho0 = self.rho0
        tau_f = self.tau_f
        
        # If rho0 is not populated then set equal ground state populations
        if(not rho0.any()):
            print("Populating ground states equally as the initial condition.")
            population = 1/len(G)
            for g in G:
                self.setRho0(g, g, population)
        
        # Resize rho_t
        self.rho_t = [ [0 for j in range(len(time))] for i in range(self.n*self.n)]
        
        if((beam_profile_averaging) and (doppler_averaging)):
            if(laser_power):
                te.timeEvolutionGaussianAndDopplerAveraging(n, E, G, Q, self.Q_decay, tau, laser_power, r_sigma, n_beam_averaging, laser_wavelength, doppler_width, doppler_detunings, time, rho0, self.rho_t, tau_f = tau_f, rabi_scaling = rabi_scaling, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else:
                print("Need to have laser_power attribute in LaserAtomSystem to use beam profile avergaing! Equate <LaserAtomSystem>.laser_power to a power in milliWatts.")
        
        elif(beam_profile_averaging):
            if(laser_power):
                te.timeEvolutionGaussianAveraging(n, E, G, Q, self.Q_decay, tau, laser_power, r_sigma, n_beam_averaging, laser_wavelength, time, rho0, self.rho_t, tau_f = tau_f, rabi_scaling = rabi_scaling, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else:
                print("Need to have laser_power attribute in LaserAtomSystem to use beam profile avergaing! Equate <LaserAtomSystem>.laser_power to the power of the laser in mW.")
            
        elif(doppler_averaging):
            if(laser_intensity):
                te.timeEvolutionDopplerAveraging(n, E, G, Q, self.Q_decay, tau, laser_intensity, laser_wavelength, doppler_width, doppler_detunings, time, rho0, self.rho_t, tau_f = tau_f, rabi_scaling = rabi_scaling, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else: 
                print("Need to have laser_intensity attribute in LaserAtomSystem! Equate <LaserAtomSystem>.laser_intensity to the intensity of the laser in mW/mm^2.")

        else:
            if(laser_intensity):
                te.timeEvolution(n, E, G, Q, self.Q_decay, tau, laser_intensity, laser_wavelength, time, rho0, self.rho_t, tau_f = tau_f, rabi_scaling = rabi_scaling, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else: 
                print("Need to have laser_intensity attribute in LaserAtomSystem! Equate <LaserAtomSystem>.laser_intensity to the intensity of the laser in mW/mm^2.")
            
            