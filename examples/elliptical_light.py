"""Elliptically polarised light on the Rb-87 D2 cycling transition (5S1/2 F=2 -> 5P3/2 F'=3).

Before LASED v2, complex polarisation amplitudes made the equations of motion non-Hermitian and the
populations grew without bound. This script checks that rho(t) stays a valid density matrix.

Run: python examples/elliptical_light.py
"""
import numpy as np
import LASED as las

wavelength = 780.241209686e-9   # m
tau = 26.2348                   # ns, 5P3/2 lifetime
G = las.generateSubStates(label_from=1, w=0, L=0, S=1/2, I=3/2, F=2)
E = las.generateSubStates(label_from=6, w=las.angularFreq(wavelength), L=1, S=1/2, J=3/2, I=3/2, F=3)

# Ellipse with major axis at 20 degrees from x and ellipticity angle 15 degrees, travelling along z
pol = las.elliptical(angle=np.deg2rad(20), ellipticity=np.deg2rad(15))
print("spherical amplitudes {q: a_q}:", {q: np.round(a, 4) for q, a in pol.items()})

rb = las.LaserAtomSystem(E, G, tau, laser_wavelength=wavelength, laser_intensity=1.0, polarisation=pol)
time = np.linspace(0, 500, 501)
rho = rb.timeEvolution(time)

hermiticity = np.abs(rho - rho.conj().transpose(0, 2, 1)).max()
trace = np.trace(rho, axis1=1, axis2=2).real
min_eig = min(np.linalg.eigvalsh(r).min() for r in rho)
print(f"max |rho - rho^dag| = {hermiticity:.1e},  trace in [{trace.min():.12f}, {trace.max():.12f}],  "
      f"min eigenvalue = {min_eig:.1e}")
print("excited fraction at 500 ns:", round(rb.excited_population[-1], 5))
print("ground populations m_F=-2..2 at 500 ns:", np.round([rb.population(g)[-1] for g in G], 4))
