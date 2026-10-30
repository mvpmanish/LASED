"""Two-photon Rydberg excitation of Rb: 5S1/2 -> 5P3/2 -> 70S1/2 (fine structure, hyperfine ignored).

Two lasers form a ladder. The intermediate state is detuned by Delta, giving an effective two-photon
Rabi frequency Omega_1 Omega_2 / (2 Delta). Spontaneous emission from 5P3/2 (scattering error) and
Rydberg decay (modelled as loss) are included, so the script prints a small error budget for a pi pulse.

Run: python examples/rydberg_two_photon.py
"""
import numpy as np
import LASED as las

twopi = 2*np.pi
lam1, lam2 = 780.241e-9, 479.6e-9                       # 5S-5P3/2 and 5P3/2-70S (approx.)
w_P = las.angularFreq(lam1)
w_R = w_P + las.angularFreq(lam2)
S = las.generateSubStates(label_from=1, w=0, L=0, S=1/2)                 # 5S1/2 m=+-1/2
P = las.generateSubStates(label_from=3, w=w_P, L=1, S=1/2, J=3/2)        # 5P3/2 m=-3/2..3/2
R = las.generateSubStates(label_from=7, w=w_R, L=0, S=1/2)               # 70S1/2 m=+-1/2

Omega1, Omega2, Delta = twopi*0.100, twopi*0.100, twopi*2.0              # Grad/s
s_up, p_top, r_up = S[1], P[3], R[1]                                     # m = +1/2, +3/2, +1/2
red = las.Laser(S, P, lam1, rabi=Omega1/abs(las.coupling(p_top, s_up, +1)),
                polarisation=las.circular(+1), detuning=Delta, name="780 nm")
blue = las.Laser(P, R, lam2, rabi=Omega2/abs(las.coupling(r_up, p_top, -1)),
                 polarisation=las.circular(-1), detuning=-Delta, name="480 nm")
system = las.AtomicSystem(S + P + R, lasers=[red, blue],
                          decays=[las.Decay(P, S, tau=26.24)],
                          losses=[las.Loss(R, tau=150e3)])            # ~150 us at 300 K

Omega_eff = Omega1*Omega2/(2*Delta)
t_pi = np.pi/Omega_eff
rho0 = np.zeros((system.n, system.n), dtype=complex)
rho0[1, 1] = 1                                                           # 5S1/2 m=+1/2
time = np.linspace(0, 2*t_pi, 401)
rho = system.evolve(time, rho0)
pops = system.populations(rho)
P_R = pops[8]
P_P = sum(pops[i] for i in range(3, 7))
k = np.argmax(P_R)
print(f"Effective two-photon Rabi frequency: {Omega_eff/twopi*1e3:.2f} MHz, pi time {t_pi:.1f} ns")
print(f"Peak Rydberg population {P_R[k]:.5f} at t = {time[k]:.1f} ns")
print(f"  population in 5P3/2 at that time:        {P_P[k]:.2e}")
print(f"  population lost (Rydberg decay + leakage): {1 - np.trace(rho[k]).real:.2e}")
print(f"  ground population in the other m state:  {pops[1][k]:.2e}  (optical pumping by scattering)")
# Error budget: photons scattered from 5P3/2 during the pi pulse = Gamma_P * integral of P_5P dt
trapezoid = getattr(np, "trapezoid", None) or np.trapz
scattered = (1/26.24)*trapezoid(P_P[:k+1], time[:k+1])
estimate = (1/26.24)*t_pi*(Omega1**2 + Omega2**2)/(8*Delta**2)
print(f"  photons scattered from 5P3/2 during the pulse: {scattered:.2e} (estimate {estimate:.2e})")
print(f"  pi-pulse infidelity 1 - P_R:                 {1 - P_R[k]:.2e}")
