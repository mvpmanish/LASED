"""
Generate Sub States.
Generates sub-states of a quantum state
"""

from LASED.state import *
import numpy as np

def generateSubStates(label_from, w, L, S, J = None, I = None, F = None):
    """Generates a vector of sub-states given quantum numbers.
    
    Parameters:
        label_from(int) : the number assigned to the first state generated - goes from -m_F to +m_F so the number given will be assigned to -m_F
        w(float) : the angular frequency (in Grad/s) given to the states generated
        L(int) : orbital angular momentum
        S(int) : spin quantum number
        J(int) : total angular momentum
        I(int) : nuclear isospin 
        F(int) : total angular momentum + isospin
    Returns:
        A vector of State objects from -m_F to +m_F with labels beginning at label_from
    """
    # Set up a vector
    states = []
    
    # Generate quantum numbers and labels
    if(J == None):
        J = L + S
    if(I == None):
        I = 0 
    if(F == None):
        F = I + J  # Estimate the value of F
        
    m_F = np.linspace(-1*F, F, 2*F+1, dtype = int)   # Generate substates
    labels = np.linspace(label_from, label_from+len(m_F)-1, len(m_F), dtype = int)  # The label_from + len(m_F) -1 gives the correct state labelling. If just used +len(m_F) would have one over the number of sub-states
    
    # Generate the states and input into a vector
    for i in range(len(m_F)):
        states.append(State(label = labels[i], w = w, m = m_F[i], L = L, S = S, J = J, I = I, F = F))
    return states
    