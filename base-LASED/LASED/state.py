'''
Class definition for an atomic state
'''

# Class to hold information to construct an atomic state
class State():
    # Parameterised constructor
    def __init__(self, label, w, m, L, S, J = None, I = None, F = None):
        self.label = label  # number labelling of state e.g. state |1> would have label 1
        self.w = w  # angular frequency of state in Grad/s
        
        self.L = L  # Orbital angular momentum quantum number
        self.S = S  # Spin quantum number
        self.m = m  # m is the degeneracy of total angular momentum
        
        if(J):
            self.J = J  # Coupling between L and S
        else:
            self.J = L+S 
        if(I):
            self.I = I  # Nuclear spin quantum number
        else:
            self.I = 0  
        if(F):
            self.F = F  # Total angular momentum, coupling between J and I
        else:
            self.F = self.J + self.I
         
        
        