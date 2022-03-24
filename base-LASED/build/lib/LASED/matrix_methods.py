"""
Defines functions to manipulate matrices.
"""

import numpy as np
import scipy.linalg as la

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
