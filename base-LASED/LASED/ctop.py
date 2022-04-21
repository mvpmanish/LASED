# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 18:58:09 2022

@author: FuckNvidia
"""

import LASED.CppLASED
from LASED.time_evolution_matrix import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j


def couplingTable_q1(E,G):
    qtable = np.empty((len(E),len(G)),dtype= np.float64)
    q = 1
    for i in range(len(E)):
        for j in range(len(G)):
            qtable[i,j] = coupling(E[i], G[j], q)

    return qtable

def couplingTable_q0(E,G):
    qtable = np.empty((len(E),len(G)),dtype= np.float64)
    q = 0
    for i in range(len(E)):
        for j in range(len(G)):
            qtable[i,j] = coupling(E[i], G[j], q)

    return qtable

def couplingTable_qn1(E,G):
    qtable = np.empty((len(E),len(G)),dtype= np.float64)
    q = -1
    for i in range(len(E)):
        for j in range(len(G)):
            qtable[i,j] = coupling(E[i], G[j], q)

    return qtable


def stateWrapping(states):
    """
    Parameters:
        states (list of state object)

    """
    n = len(states)
    statearray = np.zeros((n,8), np.float64)
    for i in range(n):
        statearray[i,0] = states[i].label
        statearray[i,1] = states[i].w
        statearray[i,2] = states[i].L
        statearray[i,3] = states[i].S
        statearray[i,4] = states[i].m
        statearray[i,5] = states[i].J
        statearray[i,6] = states[i].I
        statearray[i,7] = states[i].F

    return statearray


def wrapping(laser_atom_system, detuning = 0.0,atomic_velocity = 0.0):

    laser_intensity = laser_atom_system.laser_intensity
    laser_power = laser_atom_system.laser_power
    tau_f = laser_atom_system.tau_f
    tau_b = laser_atom_system.tau_b
    rabi_factors = laser_atom_system.rabi_factors
    rabi_scaling=laser_atom_system.rabi_scaling
    if rabi_scaling == None:
        rabi_scaling = 1
    if rabi_factors == None:
        rabi_factors = [1]
    if tau_f == None:
        tau_f = 0
    if tau_b == None:
        tau_b = 0

    if laser_intensity == None:
        laser_intensity = 0
    if laser_power == None:
        laser_power = 0

    all_states = laser_atom_system.G + laser_atom_system.E


    CppLaserAtomSystem = LASED.CppLASED.SolveLaserAtomSystem(
        stateWrapping(laser_atom_system.E),
        stateWrapping(laser_atom_system.G),
        laser_atom_system.tau,
        laser_atom_system.Q,
        laser_atom_system.laser_wavelength,
        couplingTable_qn1(all_states,all_states),
        couplingTable_q0(all_states,all_states),
        couplingTable_q1(all_states,all_states),
        float(detuning),
        float(laser_intensity),
        float(laser_power),
        float(tau_f),
        float(tau_b),
        [rabi_scaling],
        [float(i) for i in rabi_factors],
        float(atomic_velocity))
    return CppLaserAtomSystem


class wrappingsol:
    def __init__(self, sol):
        self.D = sol.get_D()
        self.V = sol.get_V()
        self.invV = sol.get_invV()

    def timeEvolution0(self, time, rho0):

        return np.matmul(self.V, np.matmul(np.diag(np.exp(self.D.transpose()[0]*time)), np.matmul(self.invV, rho0)))

    def timeEvolution1(self, timearray,rho0):
        result = np.empty((len(timearray), len(rho0)), dtype = np.complex128)
        for i in range(len(result)):
            result[i] = self.timeEvolution0(timearray[i], rho0)
        return result

    def timeEvolution(self, timearray,rho0,a,b):
        result = np.empty(len(timearray), dtype = np.complex128)
        x = a*int(np.sqrt(len(rho0)))+b

        for i in range(len(result)):
            result[i] = self.timeEvolution0(timearray[i], rho0)[x,0]
        return result

def get_S(angle, csv):
    f = np.empty(int((len(csv[0])-1)/2), dtype=np.complex128)
    for i in range(len(f)):
        f[i] = csv[angle, 2*i+1] + csv[angle, 2*i+2]*1j
    return f

def get_rho(s):
    n = len(s)
    rho0 = np.empty(n*n, dtype=np.complex128)
    Sum = 0
    for i in range(n):
        for j in range(n):
            index = j + i * n
            rho0[index] =s[i] * np.conjugate(s[j])
            if i == j:
                Sum += rho0[index]

    return np.reshape(rho0,(n, n)) / Sum

def get_rho_full(rho0, n):
    rho = np.zeros((n, n), dtype=np.complex128)
    lenth = int(np.sqrt(len(rho0)))
    #rho0 = np.reshape(rho0,(lenth,lenth))
    for i in range(lenth):
        for j in range(lenth):
            rho[i,j] = rho0[i,j]
    rho = np.reshape(rho, n*n)
    return rho


#%%
def printinC(array):
    print("{")
    for i in range(len(array)):
        line = ""
        for j in range(len(array[i])):
            line += str(array[i,j]) + ","
        print("{"+line+"},")
    print("}")


#printinC(couplingTable_q1(all_states,all_states))


def pop(sol,rho0, time0,i,j):
    n = len(time0)
    result = np.zeros(n)
    for I in range(n):
        result[I] = abs(sol.timeEvolution(rho0,time0[I])[i*int(np.sqrt(len(rho0)))+j][0])
    return result
