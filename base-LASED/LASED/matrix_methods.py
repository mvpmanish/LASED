# Function definitions for methods which operate on matrices

import numpy as np
import scipy.linalg as la

# Diagonalise a matrix.
# Input: A is a square matrix
# Output: diagonalised A
def diagonalise(A):
    eigenvals, eigenvecs = la.eig(A)
    return np.diag(eigenvals)


# Calculate a matrix of eigenvectors.
# Input: A is a square matrix
# Output: matrix of eigenvectors
def matrixOfEigenvec(A):
    eigenvals, eigenvecs = la.eig(A)
    return eigenvecs

# For debugging purposes
def printNonZeroMatrixElements(A):
    for i, row in enumerate(A):
        for j, element in enumerate(row):
            if(element != 0):
                print(f"[{i}, {j}] = {element}")