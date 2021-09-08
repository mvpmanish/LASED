"""
Defines functions to manipulate matrices.
"""

import numpy as np
import scipy.linalg as la

def diagonalise(A):
    """Diagonalise a matrix.
    
    Parameters:
        A (ndarry): A square matrix.
    
    Returns:
        (ndarry) A diagonalised matrix.
    """
    eigenvals, eigenvecs = la.eig(A)
    return np.diag(eigenvals)



def matrixOfEigenvec(A):
    """Calculate a matrix of eigenvectors.
    Parameters: 
        A (ndarray): A square matrix.
        
    Returns:
        (ndarray): A matrix of eigenvectors
    """
    eigenvals, eigenvecs = la.eig(A)
    return eigenvecs

#----------------------------------FOR DEBUGGING------------------------------
def printNonZeroMatrixElements(A):
    """Prints non-zero matrix elements
    
    Parameters:
        A (ndarray): A matrix
    """
    for i, row in enumerate(A):
        for j, element in enumerate(row):
            if(element != 0):
                print(f"[{i}, {j}] = {element}")