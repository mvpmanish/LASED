"""Benchmark: Cs-133 D2 line (Tutorial 6), 48 sub-states / 2304 coupled equations, 500 ns.

Usage: python benchmarks/bench_cs_d2.py
"""
import time
import numpy as np
import LASED as las

PI = np.pi
wl = 852.356e-9
w_e = las.angularFreq(wl)
I = 7/2
w1, w2, w3, w4 = 9.192631770*2*PI, 0.15121*2*PI, 0.20124*2*PI, 0.251*2*PI
G = las.generateSubStates(1, -w1, 0, 1/2, I=I, F=3) + las.generateSubStates(8, 0, 0, 1/2, I=I, F=4)
E = (las.generateSubStates(17, w_e-(w4+w3+w2), 1, 1/2, I=I, F=2) + las.generateSubStates(22, w_e-(w4+w3), 1, 1/2, I=I, F=3)
     + las.generateSubStates(29, w_e-w4, 1, 1/2, I=I, F=4) + las.generateSubStates(38, w_e, 1, 1/2, I=I, F=5))
t = np.linspace(0, 500, 501)


def run(label, **kw):
    s = las.LaserAtomSystem(E, G, 30.473, [0], wl, laser_intensity=50)
    tic = time.perf_counter()
    rho = s.timeEvolution(t, **kw)
    print(f"{label:<45s} {time.perf_counter() - tic:8.2f} s   excited pop(500 ns) = {s.excited_population[-1]:.6f}")
    return rho


if __name__ == "__main__":
    print(las.backendInfo())
    a = run("CW, pi light, rust (adaptive DP5)", backend="rust")
    b = run("CW, pi light, expm (exact propagator)", backend="expm")
    print("max |rust - expm| =", np.abs(a - b).max())
    s = las.LaserAtomSystem(E, G, 30.473, laser_wavelength=wl, laser_intensity=50, polarisation=las.elliptical(0.3, 0.2))
    tic = time.perf_counter(); s.timeEvolution(t, backend="rust")
    print(f"{'CW, elliptical light, rust':<45s} {time.perf_counter() - tic:8.2f} s   excited pop(500 ns) = {s.excited_population[-1]:.6f}")
    run("Gaussian pulse + chirp, rust", backend="rust", pulse=las.pulses.gaussian(t0=150, fwhm=100),
        detuning_t=las.pulses.linearChirp(2*PI*1e-3, t0=150))
    dets = np.linspace(-3, 3, 16)*2*PI*0.3
    run("Doppler average (16 classes), rust", backend="rust", doppler_averaging=True,
        doppler_width=2*PI*0.3, doppler_detunings=dets)
