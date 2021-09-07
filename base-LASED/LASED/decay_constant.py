'''
This is a file to define functions to calculate the generalised decay constants
Author: Manish Patel
Date created: 27/07/2021
'''

from LASED.state import *
from LASED.half_rabi_freq import *

def generalisedDecayConstant(ep, epp, g, G, Q_decay):
    """  Calculates the branching ratio for the generalised decay constant.
    This ratio must be multiplied by 1/tau to get the generalised decay constant in Grad/s. This is then used for evaluating vertical coherences.
    Parameters:
        ep - excited State object
        epp - excited State object
        g - ground State object
        G - list of all ground State objects
        tau - lifetime of state in ns
        Q_decay - list of decay channel polarisations allowed by transition rules, usually [-1, 0, 1]
    
    Returns:
        Value of the branching ratio for the  generalised decay constant gamma_{ep, epp, g}
    """
    # Calculate the total branching from the excited state to all ground states
    sum_decay_channels_epg = 0
    sum_decay_channels_eppg = 0
    for gp in G:
        for q in Q_decay:
            sum_decay_channels_epg += abs(coupling(ep, gp, q)*coupling(ep, gp, q))
            sum_decay_channels_eppg += abs(coupling(epp, gp, q)*coupling(epp, gp, q))

    # Calculate the decay constants separately
    gamma_epg = 0
    gamma_eppg = 0
    for q in Q_decay:
        gamma_epg += abs(coupling(ep, g, q)*coupling(ep, g, q))
        gamma_eppg += abs(coupling(epp, g, q)*coupling(epp, g, q))

    gamma_epg = gamma_epg/(sum_decay_channels_epg)
    gamma_eppg = gamma_eppg/(sum_decay_channels_eppg)

    # Calculate the generalised decay constant
    gamma_epeppg = np.power(gamma_epg*gamma_eppg, 0.5)
    # Calculate the sign: only one polarisation will result in non-zero coupling so can sum over all
    for q in Q_decay:
        if(coupling(ep, g, q)*coupling(epp, g, q) < 0): 
            gamma_epeppg = -1*gamma_epeppg
    return gamma_epeppg