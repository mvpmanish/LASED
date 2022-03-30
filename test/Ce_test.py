# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:38:56 2022

@author: FuckNvidia
"""

import LASED as las
from LASED.ctop import *
from LASED.CppLASED import *
import time
import numpy as np
import matplotlib.pyplot as plt

#%%
# 6^2S_{1/2} -> 6^2P_{3/2}
wavelength_cs = 852.356e-9  # Wavelength in nm
w_e = las.angularFreq(wavelength_cs)
tau_cs = 30.473 # in ns

I_cs = 7/2  # Isospin for sodium
PI = np.pi

# Energy Splittings
w1 = 9.192631770*2*PI # Splitting of 6^2S_{1/2}(F' = 3) -> (F' = 4) in Grad/s (Exact due to definition of the second)
w2 = 0.15121*2*PI  # Splitting between 6^2P_{3/2} F = 2 and F = 3 in Grad/s
w3 = 0.20124*2*PI  # Splitting between 6^2P_{3/2} F = 3 and F = 4 in Grad/s
w4 = 0.251*2*PI  # Splitting between 6^2P_{3/2} F = 4 and F = 5 in Grad/s

# Detunings
w_Fp3 = -1*w1
w_F2 = w_e-(w4+w3+w2)
w_F3 = w_e-(w4+w3)
w_F4 = w_e-w4
w_F5 = w_e

# Create states
# 6^2S_{1/2}
Fp3 = las.generateSubStates(label_from = 1, w = w_Fp3, L = 0, S = 1/2, I = I_cs, F = 3)
Fp4 = las.generateSubStates(label_from = 8, w = 0, L = 0, S = 1/2, I = I_cs, F = 4)

# 5^2P_{3/2}
F2 = las.generateSubStates(label_from = 17, w = w_F2, L = 1, S = 1/2, I = I_cs, F = 2)
F3 = las.generateSubStates(label_from = 22, w = w_F3, L = 1, S = 1/2, I = I_cs, F = 3)
F4 = las.generateSubStates(label_from = 29, w = w_F4, L = 1, S = 1/2, I = I_cs, F = 4)
F5 = las.generateSubStates(label_from = 38, w = w_F5, L = 1, S = 1/2, I = I_cs, F = 5)

# Declare excited and ground states
G_cs = Fp3 + Fp4
E_cs = F2 + F3 + F4 + F5

# Laser parameters
intensity_cs = 50 # mW/mm^-2
Q_cs = [0]

# Simulation parameters
start_time = 0
stop_time = 500 # in ns
time_steps = 501
time_cs = np.linspace(start_time, stop_time, time_steps)
cs_system = las.LaserAtomSystem(E_cs, G_cs, tau_cs, Q_cs, wavelength_cs, laser_intensity = intensity_cs)


population = 1/len(cs_system.G)
for g in cs_system.G:
    cs_system.setRho_0(g, g, population)

        # Resize rho_t




#%%



tic = time.perf_counter()

Csys = wrapping(cs_system, 0, 0)
sol = Csys.get_solution()


toc = time.perf_counter()
print(f"The code finished in {toc-tic:0.4f} seconds")
#%%

tic = time.perf_counter()

solution = wrappingsol(sol)


fig, ax = plt.subplots()

for i in range(38,42):
    lll=abs(solution.timeEvolution(time_cs, cs_system.rho_0,i,i))
    ax.plot(time_cs,lll,"-")


ax.set(xlabel='test', ylabel='test',
       title='test')
ax.grid()


plt.show()



toc = time.perf_counter()
print(f"The code finished in {toc-tic:0.4f} seconds")







#%%





tic = time.perf_counter()



sol.timeEvolution(cs_system.rho_0,1)

toc = time.perf_counter()

print(f"The code finished in {toc-tic:0.4f} seconds")



tic = time.perf_counter()



solution.timeEvolution0(1, cs_system.rho_0)

toc = time.perf_counter()

print(f"The code finished in {toc-tic:0.4f} seconds")

