'''
Definition of the index function
'''
from state import *


def index(i, j, n):
    '''
    Function to return the index of an element of a 2D square array of 
    dimension n and indices [i,j] if converted to a 1D array
    parameters are States and the number of States in the system
    '''
    return ((i.label-1)*n+j.label)-1  # index starts at 0


def getStateLabelsFromLineNo(line_no, n):
    '''
    Gets the indeces for labelling rho [i, j] from the line number of an array
    Input: line number of an array e.g. position of row in A or column in rho_t
    Return: A tuple of [i, j]
    '''
    division = int((line_no + 1)/n)  # Gives label i
    remainder = (line_no + 1)%n  # Gives label j
    if (remainder == 0):
        i = division
        j = n
    else:
        i = (division+1)
        j = remainder
    return (i,j)