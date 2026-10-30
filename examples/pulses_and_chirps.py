"""Time-dependent control: a Gaussian pi pulse, a Ramsey sequence with a phase step, and adiabatic
rapid passage (ARP) with a chirped detuning, on the Ca 1S0 -> 1P1 pi transition (an exact two-level system).

Run: python examples/pulses_and_chirps.py
"""
import numpy as np
import LASED as las

wavelength = 422.8e-9
G = las.generateSubStates(label_from=1, w=0, L=0, S=0)
E = las.generateSubStates(label_from=2, w=las.angularFreq(wavelength), L=1, S=0)
m0 = E[1]                                              # the m=0 excited state
c0 = abs(las.coupling(m0, G[0], 0))                   # pi-transition coupling coefficient
rabi = 2*np.pi*0.05                                    # target Rabi frequency (Grad/s = rad/ns)


def two_level(envelope=None, detuning_t=None):
    laser = las.Laser(G, E, wavelength, rabi=rabi/c0, polarisation=las.pi_polarised(),
                      envelope=envelope, detuning_t=detuning_t)
    return las.AtomicSystem(G + E, lasers=[laser])   # no decay: coherent control


# 1) Gaussian pi pulse: area = rabi * sigma * sqrt(2 pi) = pi
sigma = np.pi/(rabi*np.sqrt(2*np.pi))
pulse = las.pulses.gaussian(t0=6*sigma, sigma=sigma)
rho = two_level(pulse).evolve(np.linspace(0, 12*sigma, 3))
print(f"Gaussian pi pulse (area {pulse.area(rabi, 0, 12*sigma)/np.pi:.4f} pi): P_e = {rho[-1, 2, 2].real:.6f}")

# 2) Ramsey: pi/2 - wait - pi/2 with the second pulse's phase scanned
T = (np.pi/2)/rabi
for phase in (0, np.pi/2, np.pi):
    seq = las.pulses.square(0, T) + las.pulses.square(3*T, 4*T).withPhase(phase)
    rho = two_level(seq).evolve([0, 4*T], max_step=T/50, n_ctrl=40001)
    print(f"Ramsey, second-pulse phase {phase:.2f} rad: P_e = {rho[-1, 2, 2].real:.4f}")

# 3) ARP: sweep the detuning through resonance under a sin^2 envelope -> robust inversion
for duration in (100, 500, 3000):
    chirp = las.pulses.linearChirp(rate=6/duration, t0=duration/2)          # -3 -> +3 Grad/s
    rho = two_level(las.pulses.sinSquared(0, duration), chirp).evolve([0, duration])
    print(f"ARP over {duration:5d} ns: P_e = {rho[-1, 2, 2].real:.4f}")
