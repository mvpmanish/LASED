"""LASED: Laser-Atom interaction Simulator derived from quantum ElectroDynamics.

Quick start (see AGENTS.md / the docs for more):

    >>> import LASED as las, numpy as np
    >>> G = las.generateSubStates(label_from=1, w=0, L=0, S=0)                       # 1S0
    >>> E = las.generateSubStates(label_from=2, w=las.angularFreq(422.8e-9), L=1, S=0)  # 1P1
    >>> calcium = las.LaserAtomSystem(E, G, tau=4.6, laser_wavelength=422.8e-9, laser_intensity=100,
    ...                               polarisation=las.elliptical(angle=0, ellipticity=0.3))
    >>> rho = calcium.timeEvolution(np.linspace(0, 50, 501))    # (501, 4, 4) complex array

Main entry points:
    LaserAtomSystem   one laser driving ground states G -> excited states E (the classic LASED API)
    AtomicSystem      any number of states, lasers (Laser), decays (Decay), losses (Loss), dephasing
    pulses            pulse envelopes, phases, chirps and moving-atom trajectories
    linear, circular, elliptical, jones   polarisation helpers (spherical amplitudes {q: a_q})
    backendInfo()     whether the Rust backend is available
"""
from LASED.laser_atom_system import *
from LASED.state import *
from LASED.constants import *
from LASED.detuning import *
from LASED.half_rabi_freq import *
from LASED.decay_constant import *
from LASED.generate_sub_states import *
from LASED.index import *
from LASED.density_matrix import *
from LASED.rotation import *
from LASED.angular_shape import *
from LASED.save_to_csv import *
from LASED.time_evolution import *
from LASED.time_evolution_matrix import *
from LASED.engine import (AtomicSystem, Laser, Decay, Loss, Dephasing, Ensemble, Model,
                          dopplerEnsemble, gaussianBeamEnsemble, commutatorSuperoperator,
                          dissipatorSuperoperator)
from LASED.polarisation import (jones, linear, circular, elliptical, pi_polarised, fromRabiFactors,
                                asPolarisation)
from LASED.zeeman import landeGJ, landeGF, zeemanShift, applyZeemanShift
from LASED._backend import backendInfo, HAVE_RUST
from LASED import pulses

__version__ = "2.0.0"
