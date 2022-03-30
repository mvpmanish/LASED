# LASED: Laser Atom interaction Simulator using quantum ElectroDynamics

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
- `rotate` the laser-atom system's density matrix at t = 0, defined as `rho_0` to a different reference frame and then time evolve using the Euler angles.
- Obtain the density matrix for the excited state and ground states over all simulation time.
- Plot the time evolution of the `angularShape` of the excited or lower atomic state's electron cloud.

## Tutorials

Check out readthedocs for detailed tutorials and a guide for how to use the library: https://lased.readthedocs.io/en/latest/

## Version History

**v1.0**
- Can plot the `angularShape` of the excited or lower atomic state's electron cloud for all simulation time.
- Increased speed of the `timeEvolution` by a factor of 2.

**v0.4**:
- Ability to model decay to other states not coupled to by the laser (e.g. non-radiative decay) using the keyword `tau_b` when instantiating the laser-atom system.
- Can export the symbolic printing of the equations of motion as a .tex and/or a .pdf file file using the keyword `pretty_print_eq_tex = True` and `pretty_print_eq_pdf = True` when performing a `timeEvolution` but the keyword `pretty_print_eq_filename` must be given a string to give the new file(s) a name. Note: to export to pdf `pdflatex` must be installed on your system to convert the .tex file to a .pdf file.

## Further Reading

- A paper on this simulator can be found here: https://arxiv.org/abs/2203.12535.

## Acknowledgements

Thank you to Professor Andrew Murray, Dr Matthew Harvey, and Parinya Udommai for their continued support with this library and project.

Please cite this library if you are using it using the paper found here: https://arxiv.org/abs/2203.12535.
