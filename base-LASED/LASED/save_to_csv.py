"""
This is a file to define the function to save rho_t as a csv file.
"""

import csv
from LASED.state import *
from LASED.index import *
from itertools import zip_longest

def writeCSV(filename, headers, data, precision = None):
    """Creates a csv file using pandas dataframe.
    
    Parameters: 
        headers:  Array of strings used as column headers.
        data: Array of arrays of data used for the columns, must be same size as headers
        filename: String of the name of the csv file output
        precision (int): Precision of numbers in decimal places.
    
    Returns:
        Void.
    """
    # Check for errors
    if len(headers) != len(data):
        print("Must have the same number of headers as data columns!")
        return
    # Create csv by using zip_longest
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        rows = zip_longest(*data, fillvalue = '')
        for row in rows:
            if(precision != None):  # Set the precision
                writer.writerow([f"{i:.{precision}f}" for i in row])
            else:
                writer.writerow(row)

def stateLabel(s, state_type):
        """Returns the short-hand label of the state.
        
        Parameters:
            s (State): a State object.
            state_type (string): either 'e' or 'g' for excited state or ground state. The ground state is primed.
            
        Returns:
            "J=k;m_J=l" if k and l are J and m quantum numbers. If the state has isospin then it is "F=k;m_F=l"
        """
        
        if(s.I != None):
            letter = "F"
        else:
            letter = "J"
        
        if(state_type == "g"):  # Prime the letter if a ground state
            letter += "'"
            
        return f"{letter}={s.J};m_{letter}={s.m}"

def saveRhotAsCSV(n, E, G, time, rho_t, filename, precision = None):
    """Saves rho_t to a csv file.
    
    Creates a csv file with filename and saves each element's time evolution as a column.
    
    Parameters:
        n (int): Number of substates in the system.
        E (list of States): List of State objects in excited state of system.
        G (list of States): List of State objects in ground state of system.
        rho_t (list of lists): List of flattened 2D density matrices representing time evolution of density matrix.
        filename (string): Name of the csv file created.
        precision (int): Precision of numbers in decimal places.
    
    Returns:
        Void
    """
    
    # Create the headers: start from ground state and excited state populations then move onto coherences.
    # In the same loop sort rho_t into a format identical to the headers.
    headers = ["time(ns)"]  # First entry is time
    sorted_rho_t = []
    sorted_rho_t.append(time)  # Append the first entry
    
    g_state_type = "g"
    e_state_type = "e"
    
    for g in G:
        for gp in G:
            headers.append(f"Re({stateLabel(g, g_state_type)},{stateLabel(gp, g_state_type)})")
            sorted_rho_t.append([x.real for x in rho_t[index(g, gp, n)]])
            headers.append(f"Im({stateLabel(g, g_state_type)},{stateLabel(gp, g_state_type)})")
            sorted_rho_t.append([x.imag for x in rho_t[index(g, gp, n)]])
    for e in E:
        for ep in E:
            headers.append(f"Re({stateLabel(e, e_state_type)},{stateLabel(ep, e_state_type)})")
            sorted_rho_t.append([x.real for x in rho_t[index(e, ep, n)]])
            headers.append(f"Im({stateLabel(e, e_state_type)},{stateLabel(ep, e_state_type)})")
            sorted_rho_t.append([x.imag for x in rho_t[index(e, ep, n)]])
    for e in E:
        for g in G:
            headers.append(f"Re({stateLabel(e, e_state_type)},{stateLabel(g, g_state_type)})")
            sorted_rho_t.append([x.real for x in rho_t[index(e, g, n)]])
            headers.append(f"Im({stateLabel(e, e_state_type)},{stateLabel(g, g_state_type)})")
            sorted_rho_t.append([x.imag for x in rho_t[index(e, g, n)]])
    for g in G:
        for e in E:
            headers.append(f"Re({stateLabel(g, g_state_type)},{stateLabel(e, e_state_type)})")
            sorted_rho_t.append([x.real for x in rho_t[index(g, e, n)]])
            headers.append(f"Im({stateLabel(g, g_state_type)},{stateLabel(e, e_state_type)})")
            sorted_rho_t.append([x.imag for x in rho_t[index(g, e, n)]])


    # Write csv file
    writeCSV(filename, headers, sorted_rho_t, precision = precision)