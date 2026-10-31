"""Backend benchmark: prints a markdown table of wall-clock times for each solver.

Run on an otherwise idle machine:  python benchmarks/bench_backends.py
"""
import time
import numpy as np
import LASED as las

PI = np.pi


def cs_d2():
    """Cs-133 D2 (Tutorial 6): 48 sub-states."""
    wl, I = 852.356e-9, 7/2
    w_e = las.angularFreq(wl)
    w1, w2, w3, w4 = 9.192631770*2*PI, 0.15121*2*PI, 0.20124*2*PI, 0.251*2*PI
    G = las.generateSubStates(1, -w1, 0, 1/2, I=I, F=3) + las.generateSubStates(8, 0, 0, 1/2, I=I, F=4)
    E = (las.generateSubStates(17, w_e-(w4+w3+w2), 1, 1/2, I=I, F=2) + las.generateSubStates(22, w_e-(w4+w3), 1, 1/2, I=I, F=3)
         + las.generateSubStates(29, w_e-w4, 1, 1/2, I=I, F=4) + las.generateSubStates(38, w_e, 1, 1/2, I=I, F=5))
    return G, E, wl, 30.473


def rb87_d2():
    """Rb-87 D2: 24 sub-states."""
    wl = 780.241209686e-9
    w0, tp = las.angularFreq(wl), 2*PI
    G = (las.generateSubStates(1, -tp*6.834682611, 0, 1/2, I=3/2, F=1) + las.generateSubStates(4, 0, 0, 1/2, I=3/2, F=2))
    hf = {0: -tp*0.4958150, 1: -tp*0.4235970, 2: -tp*0.2666500, 3: 0.0}
    E, label = [], 9
    for F in range(4):
        E += las.generateSubStates(label, w0 + hf[F], 1, 1/2, J=3/2, I=3/2, F=F)
        label += 2*F + 1
    return G, E, wl, 26.2348


def timed(fn, repeat=1):
    best = np.inf
    for _ in range(repeat):
        tic = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - tic)
    return best


def main():
    print(las.backendInfo())
    rows = []
    G, E, wl, tau = cs_d2()
    t = np.linspace(0, 500, 501)

    def cs_system(pol, **laser_kw):
        laser = las.Laser(G, E, wl, intensity=50, tau=tau, polarisation=pol, **laser_kw)
        return las.AtomicSystem(G + E, [laser], [las.Decay(E, G, tau)])

    s = cs_system(las.pi_polarised())
    rows.append(("Cs D2: assemble H, jump operators and Liouvillian", "-", timed(s.build, 3)))
    for name, pol in (("pi", las.pi_polarised()), ("elliptical", las.elliptical(0.3, 0.2))):
        s = cs_system(pol)
        for b in ("expm", "rust"):
            rows.append((f"Cs D2 CW {name} light, 500 ns", b, timed(lambda: s.evolve(t, backend=b))))
    s = cs_system(las.pi_polarised())
    ens = las.dopplerEnsemble(s.lasers, sigma_v=0.25, n=16)
    for b in ("expm", "rust"):
        rows.append(("Cs D2 CW pi light, 16 Doppler classes", b, timed(lambda: s.evolve(t, ensemble=ens, backend=b))))
    s = cs_system(las.pi_polarised(), envelope=las.pulses.gaussian(t0=150, fwhm=100),
                  detuning_t=las.pulses.linearChirp(2*PI*1e-3, t0=150))
    for b in ("rust", "scipy"):
        rows.append(("Cs D2 Gaussian pulse + chirp", b, timed(lambda: s.evolve(t, backend=b))))

    G, E, wl, tau = rb87_d2()
    laser = las.Laser(G, E, wl, intensity=3, tau=tau, polarisation=las.linear(0.3),
                      envelope=las.pulses.flatTop(5, 60, rise=10).withPhase(lambda x: 0.02*x),
                      detuning_t=las.pulses.tanhChirp(2*PI*0.05, 35, 15))
    s = las.AtomicSystem(G + E, [laser], [las.Decay(E, G, tau)])
    t = np.linspace(0, 80, 9)
    for b in ("rust", "scipy"):
        rows.append(("Rb-87 D2 shaped + chirped pulse", b, timed(lambda: s.evolve(t, backend=b))))
    ens = las.dopplerEnsemble([laser], temperature=1e-3, mass=87*1.66e-27, n=16)
    for b in ("rust", "scipy"):
        rows.append(("Rb-87 D2 pulse, 16 Doppler classes", b, timed(lambda: s.evolve(t, ensemble=ens, backend=b))))

    print("\n| Problem | Backend | Time (s) |\n|---|---|---|")
    for name, b, dt in rows:
        print(f"| {name} | {b} | {dt:.3f} |")


if __name__ == "__main__":
    main()
