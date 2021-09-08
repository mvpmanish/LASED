'''
Class definition for an atomic state
'''

# Class to hold information to construct an atomic state
class State():
    """An atomic state.
    
    Attributes:
        label (int): Number labelling of state e.g. state |2> would have label 2.
        w (float): Angular frequency corresponding to energy of state in Grad/s
        L (int): Orbital angular momentum quantum number.
        S (int): Spin quantum number.
        m (int): Degeneracy of the total angular momentum.
        J (int): Resultant of L coupling to S.
        I (int): Nuclear spin quantum number.
        F (int): Total angular momentum. 
    """
    # Parameterised constructor
    def __init__(self, label, w, m, L, S, J = None, I = None, F = None):
        self.label = label  # number labelling of state e.g. state |1> would have label 1
        self.w = w  # angular frequency of state in Grad/s
        
        self.L = L  # Orbital angular momentum quantum number
        self.S = S  # Spin quantum number
        self.m = m  # m is the degeneracy of total angular momentum
        
        if(J != None):
            self.J = J  # Coupling between L and S
        else:
            self.J = L+S 
        if(I != None):
            self.I = I  # Nuclear spin quantum number
        else:
            self.I = 0  
        if(F != None):
            self.F = F  # Total angular momentum, coupling between J and I
        else:
            self.F = self.J + self.I
    
    # Print
    def __repr__(self):
        return f"State(label = {self.label}, w = {self.w}, m = {self.m}, L = {self.L}, J = {self.J}, I = {self.I}, F = {self.F})"
    
    
         
        
        