'''
The LaserAtomSystem class definition.
'''

from LASED.time_evolution import *
from LASED.density_matrix import *
from LASED.index import *
from LASED.generate_sub_states import *
from LASED.save_to_csv import *

import LASED.rotation
import numpy as np



class LaserAtomSystem:
    """A physical system composed of a laser field acting on an atomic system.
    
    Attributes: 
        Q_decay (list): List of ints describing the selection rules of the decay. Selection rules are set to +1, -1, and 0.
        rho_t (list): List of flattened 2D arrays over the time interval simulated for. This is initialised as empty as no time evolution has taken place.
        E (list): List of State objects which are the excited states of the system.
        G (list): List of State objects which are the ground states of the system.
        tau (float): Lifetime of transition in nanoseconds between excited and ground state.
        Q (list): List of laser polarisations. This can be +1, 0, -1 for right-hand circular, left-hand circular, and linear polarisation. If more than one polarisation is in the list then the system will be excited with a linear combination of the polarisations. 
        laser_wavelength (float): Wavelength of transition from ground to excited state in metres.
        laser_intensity (float): Intensity of the laser in mW/mm^2.
        laser_power (float): Power of the laser in mW. This is needed for Gaussian averaging of beam profile.
        tau_f (float): Non-radiative lifetime of the laser-atom system
        rho_0 (ndarray): 2D array creating the density matrix at t = 0.
        rabi_scaling (float): The normalisation of the Rabi frequency. The half-Rabi frequency is divided by this number. Use this if there are more than one polarisations to normalise the population.
        rabi_factors (list): Elements of this list are multiplies by the half-Rabi frequency for each polarisation. This is used to obtain elliptical polarisation or align the polarisation of the laser to a different axis e.g. rabi_factors = [1*cos(pi/2)+1*j*sin(pi/2), 2*cos(pi/2+2*j*sin(pi/2))] for Q = [1, -1] gives a more LH elliptical polarisation aligned by 45 deg to the x-axis.
    """
    
    # Class variables
    Q_decay = [1, 0 ,-1]
    rho_t = []
    time = []
    
    def __init__(self, E, G, tau, Q, laser_wavelength, laser_intensity = None, 
                 laser_power = None, tau_f = None, rabi_scaling = None, 
                rabi_factors = None):
        """
        Inits LaserAtomSystem.
        """
        self.E = E  # list of excited States
        self.G = G  # list of ground States
        self.tau = tau  # lifteime in ns/rad, N.B NIST database uses A_ki in rad/s
        self.Q = Q  # laser radiation polarisation
        self.laser_wavelength = laser_wavelength  # wavelength of the laser in nm
        self.tau_f = tau_f  # lifetime of decay to other states (can be non-radiative) in ns/rad
        self.laser_intensity = laser_intensity  # in mW/mm^2
        self.laser_power = laser_power  # in mW
        self.rho_0 = np.zeros((self.n*self.n, 1), dtype = complex)  # flattened density matrix
        self.rabi_scaling = rabi_scaling  # Normalisation factor of rabi frequencies for different polarisations
        self.rabi_factors = rabi_factors  # multiply these by rabi frequency for different polarisations
        
    @property
    def n(self):
        """ Total number of substates.
        
        Returns:
            (int) Number of substates.
        """
        return int(len(self.G)+len(self.E))
    
    @property
    def rho_e0(self):
        """ Upper state density matrix for the initial condition.
        
        Returns:
            (ndarry) Upper state density matrix composed of all states defined in E.
        """
        return getSingleStateMatrix(self.rho_0, self.n, self.E)
    
    @property
    def rho_g0(self):
        """ Lower state density matrix for the initial condition.
        
        Returns:
            (ndarry) Lower state density matrix composed of all states defined in G.
        """
        return getSingleStateMatrix(self.rho_0, self.n, self.G)
    
    @property
    def rho_et(self):
        """ Upper state density matrix for all of the time evolution.
        
        Returns:
            (list of ndarray) A list of upper state density matrices for the time evolution.
        """
        rho_et = []
        flipped_rho_t = np.transpose(self.rho_t)  # Flip to loop over all rho
        for rho in flipped_rho_t:
            new_rho = np.zeros((self.n*self.n, 1), dtype = complex)  # Placeholder
            for i, element in enumerate(new_rho):
                new_rho[i, 0] = rho[i]
            rho_et.append(getSingleStateMatrix(new_rho, self.n, self.E))
        return rho_et
    
    @property
    def rho_gt(self):
        """ Lower state density matrix for all of the time evolution.
        
        Returns:
            (list of ndarray) A list oflower state density matrices for the time evolution.
        """
        rho_et = []
        flipped_rho_t = np.transpose(self.rho_t)  # Flip to loop over all rho
        for rho in flipped_rho_t:
            new_rho = np.zeros((self.n*self.n, 1), dtype = complex)  # Placeholder
            for i, element in enumerate(new_rho):
                new_rho[i, 0] = rho[i]
            rho_et.append(getSingleStateMatrix(new_rho, self.n, self.G))
        return rho_et
    
    def __repr__(self):
        E_str = [e.label for e in self.E]
        G_str = [g.label for g in self.G]
        return f"LaserAtomSystem({E_str}, {G_str}, {self.tau}, {self.Q}, {self.Q_decay}, {self.laser_wavelength}, {self.tau_f}, {self.laser_intensity}, {self.laser_power})"
        
    
    def Rho_0(self, i, j):
        """Accessor for an element in rho_0.
        
        Parameters:
            i (State): First state index.
            j (State): Second state index.
        
        Returns:
            complex: element of the laser-atom system density matrix at t=0.
        
        Example:
            print(Rho_0(one, two))
        """
        row = index(i, j, self.n)
        return self.rho_0[row, 0]
    
    def setRho_0(self, i, j, value):
        """Sets a value to an element of rho_0.
        
        Parameters:
            i (State): First state index
            j (State): Second state index
        """
        if(value > 1):
            print("Cannot set an element of a density matrix > 1!")
            return
        else:
            row = index(i, j, self.n)
            self.rho_0[row, 0] = value
    
    def appendDensityMatrixToRho_0(self, density_rho):
        """Sets the laser-atom system density matrix at t=0 to the matrix given.
        
        Parameters:
            density_rho (ndarray): 2D array of the system density matrix.
        
        Note:
            Density matrix input must be square and the size of the matrix must match with E or G.
        """
        size = len(density_rho)
        if(size == len(G)):
            sub_states = self.G
        elif(size == len(E)):
            sub_states = self.E
        else:
            print("Size of density_rho does not match with excited or ground states")
            return
        appendDensityMatrixToFlatCoupledMatrix(self.rho_0, density_rho, sub_states, self.n)
            
    
    def clearRho_0(self):
        """Makes all values of rho_0 zero.
        """
        self.rho_0 = np.zeros((self.n*self.n, 1), dtype = complex)
    
    def Rho_t(self, i, j):
        """Accessor for an element in rho_t.
        
        Parameters:
            i (State): First state index
            j (State): Second state index
        
        Returns:
            Array of an element in laser-atom system for all of the simulation time
        
        Example:
            print(Rho_t(one, two)) prints element rho_12 if one and two are State objects corresponding 
            to label 1 and 2 respectively. 
        """
        return self.rho_t[index(i, j, self.n)]
    
    def rotateRho_0(self, alpha, beta, gamma):
        """ Rotate rho_0 by the Euler angles alpha, beta, and gamma.
        
        Parameters:
            alpha: rotation (in radians) around z-axis
            beta: rotation (in radians) about the y'-axis
            gamma: rotation (in radians) about the z''-axis 
        
        Note:
            Rotation can only be performed if isospin is zero i.e. rotation is only in J-representation
            and not in F-representation. 
        """
        # Check if any state is in F-representation
        sum_I = 0
        for e in E:
            sum_I += e.I
        for g in G:
            sum_I += g.I
        if(sum_I > 0):
            print("Warning:  Rotation can only be performed if isospin is zero i.e. rotation is only in J-representation and not in F-representation! ")
            return
        
        self.rho_0 = rotation.rotateInitialMatrix(self.rho_0, self.n, self.E, self.G, alpha, beta, gamma)
    
    def rotateRho_t(self, alpha, beta, gamma):
        """ Rotate rho_0 by the Euler angles alpha, beta, and gamma.

        Parameters:
            alpha: rotation (in radians) around z-axis
            beta: rotation (in radians) about the y'-axis
            gamma: rotation (in radians) about the z''-axis 
        
        Note:
            Rotation can only be performed if isospin is zero i.e. rotation is only in J-representation
            and not in F-representation. 
        """
        print("Optical coherences are preserved under rotation. To obtain these in a new reference frame, rotate rho_0 and then evolve in the new reference frame with the correct polarisation.")
        
        # Check if any state is in F-representation
        sum_I = 0
        for e in E:
            sum_I += e.I
        for g in G:
            sum_I += g.I
        if(sum_I > 0):
            print("Warning:  Rotation can only be performed if isospin is zero i.e. rotation is only in J-representation and not in F-representation! ")
            return
        
        rotated_rho_t = []
        # Flip to loop over all rho
        for rho in np.transpose(self.rho_t):
            new_rho = np.zeros((self.n*self.n, 1), dtype = complex)  # Placeholder
            for i, element in enumerate(new_rho):
                new_rho[i, 0] = rho[i]
            new_rho = rotation.rotateInitialMatrix(new_rho, self.n, self.E, self.G, alpha, beta, gamma)
            rotated_rho_t.append(new_rho)
        # Flip this back to the structure of rho_t
        self.rho_t = np.transpose(rotated_rho_t)[0]
                
    
    def timeEvolution(self, time, beam_profile_averaging = None, doppler_averaging = None, 
                     pretty_print_eq = None, print_eq = None, detuning = None,
                     atomic_velocity = None, r_sigma = None, n_beam_averaging = None, 
                     doppler_width = None, doppler_detunings = None):
        """ Evolves the laser-atom system over time.
        
        Produces a list of flattened 2D density matrices which correspond to the time array. This is stored in rho_t.
        
        Parameters:
            time (list): Array of the times in the time evolution in nanoseconds e.g. time = [0, 1, 2] to simulate up to 2 ns
            beam_profile_averaging (bool): Turn on averaging over a Gaussian TEM00 laser beam profile. Must have n_beam_averaging to state the number of averages to take over the beam profile and r_sigma to characterise the standard debviation of the Gaussian beam. 
            doppler_averaging (bool): Turn on averaging over the doppler profile of the atoms. Must have doppler_width and doppler_detunings as well.
            pretty_print_eq (bool): Turn on printing of the coupled differential equations symbolically using Sympy.
            print_eq (bool): Turn on printing the coupled differential equatiuons numerically.
            atomic_velocity (float): The velocity component of the atoms in the direction of the laser beam in metres/second. This is used for doppler shifting the energy levels.
            r_sigma (float): The radial distance to the 2D standard deviation in millimetres of the Gaussian beam profile of the laser.
            n_beam_averaging (int): The number of times the beam profile will be split to average over. The higher the number, the more accurate the beam profile averaging but it will be slower.
            doppler_width (float): The doppler width of the beam profile in Grad/s.
            doppler_detunings (list): List of the doppler detunings creating a doppler profile. This list will be averaged over. Must go from a -ve detuning to a +ve detuning. All detunings are in Grad/s.
            
            Note:
                Must have laser_power attribute for Gaussian averaging of the beam and must have laser_intensity attribute when not Gaussian averaging. Also, cannot print equations in doppler or Gaussian averaging. 
        """
        self.time = time  # Save the time
        
        n = self.n
        E = self.E
        G = self.G
        Q = self.Q
        tau = self.tau
        laser_power = self.laser_power
        laser_intensity = self.laser_intensity
        laser_wavelength = self.laser_wavelength
        rho_0 = self.rho_0
        tau_f = self.tau_f
        rabi_scaling = self.rabi_scaling
        rabi_factors = self.rabi_factors
        
        # If rho0 is not populated then set equal ground state populations
        if(not rho_0.any()):
            print("Populating ground states equally as the initial condition.")
            population = 1/len(G)
            for g in G:
                self.setRho_0(g, g, population)
        
        # Resize rho_t
        self.rho_t = [ [0 for j in range(len(time))] for i in range(self.n*self.n)]
        
        if((beam_profile_averaging) and (doppler_averaging)):
            if(laser_power):
                timeEvolutionGaussianAndDopplerAveraging(n, E, G, Q, self.Q_decay, tau, laser_power, r_sigma, n_beam_averaging, laser_wavelength, doppler_width, doppler_detunings, time, rho_0, self.rho_t, tau_f = tau_f, detuning = detuning, rabi_scaling = rabi_scaling, rabi_factors = rabi_factors, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else:
                print("Need to have laser_power attribute in LaserAtomSystem to use beam profile avergaing! Equate <LaserAtomSystem>.laser_power to a power in milliWatts.")
        
        elif(beam_profile_averaging):
            if(laser_power):
                timeEvolutionGaussianAveraging(n, E, G, Q, self.Q_decay, tau, laser_power, r_sigma, n_beam_averaging, laser_wavelength, time, rho_0, self.rho_t, tau_f = tau_f, detuning = detuning, rabi_scaling = rabi_scaling, rabi_factors = rabi_factors, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else:
                print("Need to have laser_power attribute in LaserAtomSystem to use beam profile avergaing! Equate <LaserAtomSystem>.laser_power to the power of the laser in mW.")
            
        elif(doppler_averaging):
            if(laser_intensity):
                timeEvolutionDopplerAveraging(n, E, G, Q, self.Q_decay, tau, laser_intensity, laser_wavelength, doppler_width, doppler_detunings, time, rho_0, self.rho_t, tau_f = tau_f, detuning = detuning, rabi_scaling = rabi_scaling, rabi_factors = rabi_factors, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else: 
                print("Need to have laser_intensity attribute in LaserAtomSystem! Equate <LaserAtomSystem>.laser_intensity to the intensity of the laser in mW/mm^2.")

        else:
            if(laser_intensity):
                timeEvolution(n, E, G, Q, self.Q_decay, tau, laser_intensity, laser_wavelength, time, rho_0, self.rho_t, tau_f = tau_f, detuning = detuning, rabi_scaling = rabi_scaling, rabi_factors = rabi_factors, print_eq = print_eq, pretty_print_eq = pretty_print_eq, atomic_velocity = atomic_velocity)
            else: 
                print("Need to have laser_intensity attribute in LaserAtomSystem! Equate <LaserAtomSystem>.laser_intensity to the intensity of the laser in mW/mm^2.")
                

    def saveToCSV(self, filename, precision = None):
        """Saves rho_t as a csv file
        
        Parameters:
            filename (string): Name of the csv file created.
            precision (int): Precision of numbers in decimal places.
        
        Returns:
            Void.
        """
        saveRhotAsCSV(self.n, self.E, self.G, self.time, self.rho_t, filename, precision = None)
        
            