# LASED: Laser Atom interaction Simulator using quantum ElectroDynamics

THIS PROJECT IS STILL IN BETA. PLEASE GIVE SOME FEEDBACK TO github.com/mvpmanish.

Many experiments using atoms and lasers are performed in physics which require knowledge and modelling about the excited state of the atomic species being studied. Steady-state models can be used to get the final equilibrium of a laser-atom system but a large number of laser-atom interactions are short-lived an decay quickly. Most models using the Louiville equation to capture the dynamics of the interaction do not use a full quantum electrodynamic picture to evolve the system over time but instead use a semi-classical approach. In this simulator all dynamics are calculated by deriving the equations from field operators. This gives a more physcially accurate model.

## Installation

Run the following to install:
```
pip install LASED
```

The source code can be found at https://github.com/mvpmanish/LASED.

## Usage

In this simulator a user defines a `State` object with all quantum numbers defined. The user then creates two vectors: one containing all the ground states and one for the excited states. The user can then define a `LaserAtomSystem` object with a laser power (or intensity) and the laser wavelength. With this object the user can:
- Use `timeEvolution` to time evolve the laser-atom system and access the time evolution of the density matrix elements over time using `Rho_t`. Can simulate very simple systems such as magnesium and calcium with no hyperfine structure to atoms with hyperfine structure and a large number of states such as caesium.
- `rotate` the laser-atom system's density matrix at t = 0, defined as `rho_0` to a different reference frame and then time evolve using the Euler angles
- Obtain the density matrix for the excited state and ground states over all simulation time

## Tutorials

Check out readthedocs for detailed tutorials and a guide for how to use the library. 

## Acknowledgements

Thank you to Professor Andrew Murray, Dr Matthew Harvey, and Parinya Udommai for their continued support with this library and project.

Please cite this library if you are using it.
