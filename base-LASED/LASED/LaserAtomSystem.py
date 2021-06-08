'''
The LaserAtomSystem class definition.
'''

from time_evolution import *

class LaserAtomSystem:
    """A user-defined laser field acting on a user-defined atomic system
    """
    
    # Class variables
    Q_decay = [1, 0 ,-1]
    
    def __init__(self, E, G, tau, Q, laser_wavelength, time, tau_f = None):
        self.E = E
        self.G = G
        self.n = int(len(G)+ len(E))
        self.laser_wavelength = laser_wavelength
        self.tau_f = tau_f