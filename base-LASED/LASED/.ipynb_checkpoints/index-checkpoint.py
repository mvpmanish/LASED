'''
Definition of the index function
'''
from LASED.state import *


def index(i, j, n):
    """Function to return the index of an element of a 2D density matrix if flattened to a 1D array
    
    Parameters:
        i (State): To get the row number of the element to be selected
        j (State): To get the column number of the element to be selected
        n (int): dimension of the square array i.e. number of states in the system
    
    Returns:
        int: Index of the element in a 1D array
    
    Example:
        To get the index of rho_44 in a 4-level system call index(four, four, 4) where "four" is the variable which stores the State object labelled 4.
    """ 
    
    return ((i.label-1)*n+j.label)-1  # index starts at 0


def getStateLabelsFromLineNo(line_no, n):
    """ Gets the indeces for labelling rho [i, j] from the line number of an array.
    
    Parameters:
        line_no (int): line number of an array e.g. position of a row in a matrix or column in rho_t
        n: the number of states in the system
    
    Returns:
        tuple: A tuple of the indeces for the matrix position
    """
    
    division = int((line_no + 1)/n)  # Gives label i
    remainder = (line_no + 1)%n  # Gives label j
    if (remainder == 0):
        i = division
        j = n
    else:
        i = (division+1)
        j = remainder
    return (i,j)