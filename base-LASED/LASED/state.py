'''
Class definition for an atomic state
'''

# Class to hold information to construct an atomic state
class State():
    # Parameterised constructor
    def __init__(self, L, S, I, m, w, label):
        self.L = L  # Orbital angular momentum quantum number
        self.S = S  # Spin quantum number
        self.m = m  # m is the degeneracy of total angular momentum
        self.J = L+S  # Total angular momentum
        self.I = I  # Nuclear spin quantum number
        self.F = I + self.J
        self.w = w  # angular frequency of state
        self.label = label  # number labelling of state e.g. state |1> would have label 1