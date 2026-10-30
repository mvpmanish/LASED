'''
This file contains the function to calculate the time evolution of the density matrix
for an atomic system interacting with a laser.

These functions keep the original LASED (v1) call signatures and are thin wrappers around the
general engine in ``LASED.engine``. New code should use ``LaserAtomSystem.timeEvolution`` or
``AtomicSystem.evolve`` directly.
'''

from LASED.constants import *
from LASED.half_rabi_freq import gaussianIntensity
from LASED.engine import AtomicSystem, Laser, Decay, Loss, Ensemble
from LASED.polarisation import fromRabiFactors

import numpy as np


def buildLegacySystem(E, G, Q, tau, laser_wavelength, laser_intensity, tau_f=None, tau_b=None,
                      detuning=None, rabi_scaling=None, rabi_factors=None, polarisation=None,
                      envelope=None, detuning_t=None, linewidth=0.0):
    """Builds an ``AtomicSystem`` for one laser driving G -> E with spontaneous decay E -> G.

    Returns:
        AtomicSystem
    """
    pol = polarisation if polarisation is not None else fromRabiFactors(Q, rabi_factors, rabi_scaling)
    laser = Laser(lower=G, upper=E, wavelength=laser_wavelength, intensity=laser_intensity, tau=tau,
                  polarisation=pol, detuning=0.0 if detuning is None else detuning,
                  envelope=envelope, detuning_t=detuning_t, linewidth=linewidth)
    losses = []
    if tau_f is not None:
        losses.append(Loss(E, tau_f))
    if tau_b is not None:
        losses.append(Loss(G, tau_b))
    return AtomicSystem(list(G) + list(E), lasers=[laser], decays=[Decay(E, G, tau)], losses=losses)


def legacyDopplerEnsemble(doppler_width, doppler_detunings):
    """Doppler ensemble with the weighting used by LASED v1 (Gaussian of width ``doppler_width`` in
    Grad/s sampled at ``doppler_detunings``)."""
    dets = np.asarray(doppler_detunings, dtype=float)
    scaler = 1/np.sqrt(2*PI*doppler_width**2)
    spacing = np.abs(np.diff(dets))
    d_doppler = np.concatenate(([spacing[0]], spacing)) if len(dets) > 1 else np.ones(1)
    weights = d_doppler*scaler*np.exp(-0.5*(dets/doppler_width)**2)
    return Ensemble(weights, dets[:, None], np.ones((len(dets), 1)))


def legacyBeamEnsemble(laser_power, r_sigma, n_intensity):
    """Gaussian beam-profile ensemble as in LASED v1: rings at radii linspace(r_sigma, 3 r_sigma, n)
    with weights (2k+1)/n^2. Returns (peak intensity in mW/mm^2, Ensemble)."""
    R = np.linspace(r_sigma, 3*r_sigma, n_intensity)
    I_peak = gaussianIntensity(laser_power, r_sigma, 0)
    I = gaussianIntensity(laser_power, r_sigma, R)
    weights = (2*np.arange(n_intensity) + 1)/n_intensity**2
    return I_peak, Ensemble(weights, np.zeros((n_intensity, 1)), np.sqrt(I/I_peak)[:, None])


def _store(rho_output, rho_t):
    """Copies an (T, n, n) result into the legacy (n*n) x T output container."""
    flat = rho_t.reshape(rho_t.shape[0], -1).T
    for i in range(flat.shape[0]):
        rho_output[i][:] = flat[i]


def _rho0_matrix(rho0, n):
    return np.asarray(rho0, dtype=complex).reshape(n, n)


def timeEvolution(n, E, G, Q, Q_decay, tau, laser_intensity, laser_wavelength, time, rho0, rho_output,
                    tau_f = None, tau_b = None, detuning = None, rabi_scaling = None, rabi_factors = None,
                    print_eq = None, atomic_velocity = None, pretty_print_eq = None, pretty_print_eq_tex = None,
                    pretty_print_eq_pdf = None, pretty_print_eq_filename = None, **kwargs):
    """Calculates the time evolution of a laser-atom system.

    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output.
    """
    system = buildLegacySystem(E, G, Q, tau, laser_wavelength, laser_intensity, tau_f, tau_b, detuning,
                               rabi_scaling, rabi_factors)
    rho_t = system.evolve(time, _rho0_matrix(rho0, n), velocity=atomic_velocity or 0.0, **kwargs)
    _store(rho_output, rho_t)
    return rho_t


def timeEvolutionDopplerAveraging(n, E, G, Q, Q_decay, tau, laser_intensity, laser_wavelength, doppler_width,
                                    doppler_detunings, time, rho0, rho_output, tau_f = None, tau_b = None,
                                    detuning = None, rabi_scaling = None, rabi_factors = None, print_eq = None,
                                    atomic_velocity = None,  pretty_print_eq = None, pretty_print_eq_tex = None,
                                    pretty_print_eq_pdf = None, pretty_print_eq_filename = None, **kwargs):
    """Calculates the time evolution of a laser-atom system with a Gaussian doppler profile for the atoms.

    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output.
    """
    system = buildLegacySystem(E, G, Q, tau, laser_wavelength, laser_intensity, tau_f, tau_b, detuning,
                               rabi_scaling, rabi_factors)
    ensemble = legacyDopplerEnsemble(doppler_width, doppler_detunings)
    rho_t = system.evolve(time, _rho0_matrix(rho0, n), ensemble=ensemble, velocity=atomic_velocity or 0.0, **kwargs)
    _store(rho_output, rho_t)
    return rho_t


def timeEvolutionGaussianAveraging(n, E, G, Q, Q_decay, tau, laser_power, r_sigma, n_intensity, laser_wavelength,
                                    time, rho0, rho_output, tau_f = None, tau_b = None, detuning = None, rabi_scaling = None,
                                    rabi_factors = None, print_eq = None,  atomic_velocity = None,  pretty_print_eq = None,
                                    pretty_print_eq_tex = None, pretty_print_eq_pdf = None, pretty_print_eq_filename = None,
                                    **kwargs):
    """Calculates the time evolution of a laser-atom system with a Gaussian laser beam profile.

    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output.
    """
    I_peak, ensemble = legacyBeamEnsemble(laser_power, r_sigma, n_intensity)
    system = buildLegacySystem(E, G, Q, tau, laser_wavelength, I_peak, tau_f, tau_b, detuning,
                               rabi_scaling, rabi_factors)
    rho_t = system.evolve(time, _rho0_matrix(rho0, n), ensemble=ensemble, velocity=atomic_velocity or 0.0, **kwargs)
    _store(rho_output, rho_t)
    return rho_t


def timeEvolutionGaussianAndDopplerAveraging(n, E, G, Q, Q_decay, tau, laser_power, r_sigma, n_intensity, laser_wavelength,
                                            doppler_width, doppler_detunings, time, rho0, rho_output, tau_f = None,
                                            tau_b = None, detuning = None, rabi_scaling = None, rabi_factors = None,
                                            print_eq = None, atomic_velocity = None,  pretty_print_eq = None,
                                            pretty_print_eq_tex = None, pretty_print_eq_pdf = None,
                                            pretty_print_eq_filename = None, **kwargs):
    """Calculates the time evolution of a laser-atom system with a Gaussian doppler profile for the atoms and a Gaussian laser beam profile.

    Uses a flattened density matrix rho0 and calculates the time evolution over the time specified.
    The density matrix at each time step is stored in rho_output.
    """
    I_peak, beam = legacyBeamEnsemble(laser_power, r_sigma, n_intensity)
    ensemble = legacyDopplerEnsemble(doppler_width, doppler_detunings)*beam
    system = buildLegacySystem(E, G, Q, tau, laser_wavelength, I_peak, tau_f, tau_b, detuning,
                               rabi_scaling, rabi_factors)
    rho_t = system.evolve(time, _rho0_matrix(rho0, n), ensemble=ensemble, velocity=atomic_velocity or 0.0, **kwargs)
    _store(rho_output, rho_t)
    return rho_t
