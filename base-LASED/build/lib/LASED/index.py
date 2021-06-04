'''
Definition of the index function
'''
from state import *

# Function to return the index of an element of a 2D square array of 
# dimension n and indices [i,j] if converted to a 1D array
# parameters are States and the number of States in the system
def index(i, j, n):
    return ((i.label-1)*n+j.label)-1  # index starts at 0