Pulses, multiple lasers and ensembles
=====================================

Pulses and chirps
-----------------

Any laser can carry a complex envelope :math:`A(t)` (amplitude and phase) and a time-dependent detuning
:math:`\delta(t)`::

    import LASED as las, numpy as np
    pulse = las.pulses.gaussian(t0=50, fwhm=10).withPhase(np.pi/2)
    chirp = las.pulses.linearChirp(rate=2*np.pi*1e-3, t0=50)        # Grad/s per ns
    system.timeEvolution(time, pulse=pulse, detuning_t=chirp)       # LaserAtomSystem
    # or Laser(..., envelope=pulse, detuning_t=chirp) in an AtomicSystem

Available shapes are ``constant``, ``square``, ``gaussian``, ``sinSquared``, ``blackman``,
``flatTop`` and ``fromSamples``. Pulses combine with ``*`` and ``+`` and support ``withPhase``,
``delayed`` and ``area``. Detuning helpers are ``linearChirp``, ``tanhChirp`` and
``phaseFromDetuning``. Envelopes are tabulated on ``n_ctrl`` points and linearly interpolated.

Several lasers
--------------

:class:`LASED.engine.AtomicSystem` accepts any number of :class:`~LASED.engine.Laser` objects. Each one
couples a list of lower states to a list of upper states, which covers Λ, V, ladder (Rydberg
two-photon) and tripod schemes. Decay channels (:class:`~LASED.engine.Decay`), losses
(:class:`~LASED.engine.Loss`) and dephasing (:class:`~LASED.engine.Dephasing`) are listed separately.
See ``examples/rydberg_two_photon.py``.

Ensembles and moving atoms
--------------------------

:class:`~LASED.engine.Ensemble` holds weights, per-laser detuning offsets and per-laser Rabi factors.
``dopplerEnsemble`` builds a thermal velocity distribution and ``gaussianBeamEnsemble`` a beam profile.
Ensembles combine with ``*`` (outer product). Members are integrated in parallel by the Rust backend.

``pulses.movingAtom(laser, trajectory, waist)`` converts a classical trajectory into the envelope
(position in a Gaussian beam) and the Doppler detuning that the atom sees.
