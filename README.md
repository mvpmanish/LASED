# LASED: Laser Atom interaction Simulator using quantum ElectroDynamics


Many experiments using atoms and lasers are performed in physics which require knowledge and modelling about the excited state of the atomic species being studied. Steady-state models can be used to get the final equilibrium of the laser-atom system but a large number of laser-atom interactions are short-lived an decay quickly. Most models using the Louiville equation to capture the dynamics of the interaction do not use a full quantum electrodynamic picture to evolve the system over time but instead use a semi-classical approach. In this simulator all dynamics are calculated by deriving the equations from field operators. This gives a more physcially accurate model. 

## Installation

Open up your terminal and ruun the following to install:
```python
pip install LASED
```

The source code can be found at https://github.com/mvpmanish/LASED.

## Usage

In this simulator a user defines a `State` object with all quantum numbers defined. The user then creates two vectors: one containing all the ground states and one for the excited states. The user can then define a `LaserAtomSystem` object with a laser power (or intensity) and the laser wavelength. With this object the user can:
- `timeEvolve` the laser-atom system and access the density matrix elements evolving over time as `rho_t`
- `rotate` the laser-atom system's density matrix at t = 0, defined as `rho_0` to a different reference frame and then time evolve using the Euler angles
- Obtain the density matrix for the excited state and ground states over all simulation time

The Jupyter notebooks contain usage cases of the library for various atoms.

Please cite this library if you are using it. 
