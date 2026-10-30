"""Ensembles and motion: thermal Doppler averaging of a Rabi oscillation, and an atom flying through a
focused Gaussian beam (time-dependent Rabi frequency and Doppler shift from a trajectory).

Run: python examples/doppler_and_moving_atoms.py
"""
import numpy as np
import LASED as las

wavelength = 422.8e-9
G = las.generateSubStates(label_from=1, w=0, L=0, S=0)
E = las.generateSubStates(label_from=2, w=las.angularFreq(wavelength), L=1, S=0)
mass_ca = 40*1.6605e-27

# Thermal dephasing of Rabi oscillations (1D Maxwell-Boltzmann along the beam), run in parallel in Rust
laser = las.Laser(G, E, wavelength, rabi=2*np.pi*0.2, polarisation=las.pi_polarised())
system = las.AtomicSystem(G + E, lasers=[laser], decays=[las.Decay(E, G, tau=4.6)])
time = np.linspace(0, 30, 7)
for T in (1e-3, 1.0, 300.0):          # kelvin
    ens = las.dopplerEnsemble([laser], temperature=T, mass=mass_ca, n=61)
    rho = system.evolve(time, ensemble=ens)
    print(f"T = {T:8.0e} K: P_e(t) = {np.round(rho[:, 2, 2].real, 3)}")

# An atom crossing a 20 um waist beam at 2 m/s at 60 degrees to the beam axis: it sees a Gaussian
# intensity pulse and a (constant) Doppler shift. Any trajectory works, e.g. an accelerating one.
v = 2.0
direction = np.array([np.cos(np.deg2rad(60)), 0, np.sin(np.deg2rad(60))])
trajectory = lambda t: np.outer((t - 15000)*1e-9*v, direction)           # t in ns -> metres
env, doppler = las.pulses.movingAtom(laser, trajectory, waist=20e-6)
moving = las.Laser(G, E, wavelength, rabi=2*np.pi*0.01, polarisation=las.pi_polarised(),
                   envelope=env, detuning_t=doppler)
flyby = las.AtomicSystem(G + E, lasers=[moving], decays=[las.Decay(E, G, tau=4.6)])
t = np.linspace(0, 30000, 7)
rho = flyby.evolve(t, max_step=10.0)
print("Doppler shift seen during flight (MHz):", np.round(doppler(t)/(2*np.pi)*1e3, 3))
print("Flyby: P_e(t) =", np.round(rho[:, 2, 2].real, 5))
