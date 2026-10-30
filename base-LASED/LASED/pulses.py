'''
Pulse envelopes, phases and chirps for time-dependent lasers.

A ``Pulse`` is a vectorised function of time t (ns) returning a complex amplitude A(t) that multiplies
a laser's Rabi frequency. Pulses compose: ``p1 * p2``, ``p1 + p2``, ``2.0 * p`` and ``p.withPhase(phi)``.
Detuning functions (chirps) return real angular frequencies in Grad/s and are passed to
``Laser(detuning_t=...)`` or ``LaserAtomSystem.timeEvolution(detuning_t=...)``.

Example:
    >>> import LASED as las, numpy as np
    >>> pi_pulse = las.pulses.gaussian(t0=50, fwhm=10)                 # peak amplitude 1 at t = 50 ns
    >>> shaped = las.pulses.blackman(t_start=0, duration=200).withPhase(np.pi/2)
    >>> arp = las.pulses.linearChirp(rate=2*np.pi*0.01, t0=100)         # 10 MHz/ns sweep through resonance
    >>> laser = las.Laser(lower=G, upper=E, wavelength=780e-9, rabi=2*np.pi*0.05,
    ...                   envelope=pi_pulse, detuning_t=arp)
'''

import numpy as np


class Pulse:
    """A complex envelope A(t). Call it with an array of times in ns."""

    def __init__(self, fn, description="pulse"):
        self._fn = fn
        self.description = description

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        return np.asarray(self._fn(t), dtype=complex)*np.ones_like(t, dtype=complex)

    def __mul__(self, other):
        if isinstance(other, Pulse):
            return Pulse(lambda t: self(t)*other(t), f"({self.description})*({other.description})")
        return Pulse(lambda t: other*self(t), f"{other}*({self.description})")

    __rmul__ = __mul__

    def __add__(self, other):
        if isinstance(other, Pulse):
            return Pulse(lambda t: self(t) + other(t), f"({self.description})+({other.description})")
        return Pulse(lambda t: self(t) + other, f"({self.description})+{other}")

    __radd__ = __add__

    def withPhase(self, phase):
        """Multiplies by exp(i*phase); ``phase`` is a number (rad) or a function of t (ns) returning rad."""
        if callable(phase):
            return Pulse(lambda t: self(t)*np.exp(1j*np.asarray(phase(t), dtype=float)),
                         f"({self.description})*exp(i*phase(t))")
        return Pulse(lambda t: self(t)*np.exp(1j*phase), f"({self.description})*exp(i*{phase})")

    def delayed(self, dt):
        """The same pulse shifted later in time by ``dt`` ns."""
        return Pulse(lambda t: self(np.asarray(t) - dt), f"({self.description}) delayed {dt} ns")

    def area(self, rabi, t_start, t_stop, n=20001):
        """Pulse area theta = integral of rabi*|A(t)| dt (radians) between t_start and t_stop (ns)."""
        t = np.linspace(t_start, t_stop, n)
        return float(np.trapezoid(rabi*np.abs(self(t)), t)) if hasattr(np, "trapezoid") else \
            float(np.trapz(rabi*np.abs(self(t)), t))

    def __repr__(self):
        return f"Pulse({self.description})"


def constant(amplitude=1.0):
    """A constant (continuous-wave) amplitude."""
    return Pulse(lambda t: amplitude*np.ones_like(t), f"constant({amplitude})")


def square(t_on, t_off, amplitude=1.0):
    """Square pulse: ``amplitude`` for t_on <= t < t_off, else 0. (Tabulated envelopes are linearly
    interpolated, so the edges have a rise time of about one control-grid spacing.)"""
    return Pulse(lambda t: np.where((t >= t_on) & (t < t_off), amplitude, 0.0), f"square({t_on}, {t_off})")


def gaussian(t0, fwhm=None, sigma=None, amplitude=1.0):
    """Gaussian amplitude envelope centred on t0, with full width at half maximum ``fwhm`` of the
    amplitude (or standard deviation ``sigma``) in ns."""
    if sigma is None:
        if fwhm is None:
            raise ValueError("Give fwhm or sigma")
        sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    return Pulse(lambda t: amplitude*np.exp(-0.5*((t - t0)/sigma)**2), f"gaussian(t0={t0}, sigma={sigma})")


def sinSquared(t_start, duration, amplitude=1.0):
    """sin^2 envelope from t_start to t_start + duration."""
    def fn(t):
        x = (t - t_start)/duration
        return np.where((x >= 0) & (x <= 1), amplitude*np.sin(np.pi*x)**2, 0.0)
    return Pulse(fn, f"sin^2(t_start={t_start}, duration={duration})")


def blackman(t_start, duration, amplitude=1.0):
    """Blackman window envelope (low spectral leakage) from t_start to t_start + duration."""
    def fn(t):
        x = (t - t_start)/duration
        w = 0.42 - 0.5*np.cos(2*np.pi*x) + 0.08*np.cos(4*np.pi*x)
        return np.where((x >= 0) & (x <= 1), amplitude*w, 0.0)
    return Pulse(fn, f"blackman(t_start={t_start}, duration={duration})")


def flatTop(t_start, duration, rise, amplitude=1.0):
    """Flat-top pulse with sin^2 rise and fall of length ``rise`` ns; total length ``duration``."""
    if 2*rise > duration:
        raise ValueError("2*rise must not exceed duration")

    def fn(t):
        x = t - t_start
        up = np.sin(0.5*np.pi*np.clip(x/rise, 0, 1))**2 if rise > 0 else (x >= 0).astype(float)
        down = np.sin(0.5*np.pi*np.clip((duration - x)/rise, 0, 1))**2 if rise > 0 else (x <= duration).astype(float)
        return np.where((x >= 0) & (x <= duration), amplitude*np.minimum(up, down), 0.0)
    return Pulse(fn, f"flatTop(t_start={t_start}, duration={duration}, rise={rise})")


def fromSamples(t_samples, values):
    """Envelope linearly interpolated from samples (complex values allowed), zero outside the range."""
    ts = np.asarray(t_samples, dtype=float)
    vs = np.asarray(values, dtype=complex)

    def fn(t):
        re = np.interp(t, ts, vs.real, left=0.0, right=0.0)
        im = np.interp(t, ts, vs.imag, left=0.0, right=0.0)
        return re + 1j*im
    return Pulse(fn, f"samples({len(ts)})")


# ------------------------------------------------------------------ detuning functions (Grad/s)

def linearChirp(rate, t0=0.0, offset=0.0):
    """Detuning delta(t) = offset + rate*(t - t0) in Grad/s (rate in Grad/s per ns)."""
    return lambda t: offset + rate*(np.asarray(t, dtype=float) - t0)


def tanhChirp(amplitude, t0, width):
    """Detuning delta(t) = amplitude*tanh((t - t0)/width) in Grad/s (used for robust adiabatic passage)."""
    return lambda t: amplitude*np.tanh((np.asarray(t, dtype=float) - t0)/width)


def phaseFromDetuning(detuning, t_ref=0.0, n=20001, t_max=None):
    """Converts a detuning function delta(t) into the equivalent laser phase phi(t) = -integral delta dt.

    Useful when a phase-modulated pulse (e.g. a time-optimal gate) is specified as a phase profile.
    The integral is tabulated on [t_ref, t_max] and interpolated."""
    if t_max is None:
        raise ValueError("Give t_max")
    ts = np.linspace(t_ref, t_max, n)
    d = np.asarray(detuning(ts), dtype=float)
    phi = -np.concatenate(([0.0], np.cumsum(0.5*(d[1:] + d[:-1])*np.diff(ts))))
    return lambda t: np.interp(t, ts, phi)


# ------------------------------------------------------------------ moving atoms

def movingAtom(laser, trajectory, waist, beam_centre=(0.0, 0.0, 0.0)):
    """Envelope and detuning seen by an atom moving through a Gaussian laser beam.

    Parameters:
        laser (Laser): The laser (its direction and wavelength set the Doppler shift).
        trajectory (callable): t (ns, array) -> positions (array of shape (len(t), 3)) in metres.
        waist (float): 1/e^2 intensity radius of the beam in metres.
        beam_centre (3-vector): A point on the beam axis in metres.

    Returns:
        (Pulse, callable): The Rabi envelope exp(-r_perp^2/waist^2) and the Doppler detuning
        delta(t) = -k.v(t) (Grad/s), with v obtained by differentiating the trajectory. Pass them as
        ``envelope`` and ``detuning_t`` (add any other chirp to the detuning yourself).
    """
    k_hat = np.asarray(laser.direction, dtype=float)
    k_hat = k_hat/np.linalg.norm(k_hat)
    c = np.asarray(beam_centre, dtype=float)

    def positions(t):
        return np.atleast_2d(np.asarray(trajectory(np.atleast_1d(np.asarray(t, dtype=float))), dtype=float))

    def envelope(t):
        r = positions(t) - c
        r_perp = r - np.outer(r @ k_hat, k_hat)
        return np.exp(-np.sum(r_perp**2, axis=1)/waist**2)

    def detuning(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        h = 1e-3  # ns, central difference
        v = (positions(t + h) - positions(t - h))/(2*h*1e-9)  # m/s
        return -laser.wavenumber*(v @ k_hat)*1e-9
    return Pulse(envelope, "movingAtom envelope"), detuning
