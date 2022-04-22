================
Getting Started
================

Getting Started
================

What is LASED?
---------------

LASED stands for Laser-Atom interaction Simulator derived from quantum ElectroDynamics. LASED is a python library which can:

* Calculate the time evolution of an atomic system interacting with a laser.
* Generate the equations of motion of an atom-laser system.
* Rotate an atomic system to a different reference frame.
* Calculate the time evolution of the angular shape of an atomic state.

LASED can simulate any atomic system. The sub-states, angular momenta, spin, and energies of the system are provided by the user to simulate the atomic system.

Laser parameters also need to be specified to build the laser-atom system. These parameters include detuning from the transition frequency, polarisation, and laser intensity/power.

LASED can simulate a Gaussian beam profile and Doppler averaging over the atoms to provide a more accurate model of the atom-laser system. LASED can also simulate the angular shape of an atomic state over time as it is excited by a laser. This is useful in many experiments using atoms and lasers.

Installation
--------------

You can easily install LASED by opening up a terminal and running::

  pip install LASED

The source code can be viewed `here <https://github.com/mvpmanish/LASED>`__

Using LASED
------------

Start by going to :doc:`tutorials`.

The first tutorial is a guide on how to simulate one of the simplest excitations of an atom with a laser and shows how to use most of the functionality of LASED.

The second tutorial shows how to simulate decays to other states not in the laser-excitation manifold. It also introduces rotating the system to different reference frames using the Wigner rotation matrix and hence introduces simulating a polarisation angle.

The third tutorial simulates a metastable Helium atom with some real data on what the initial state looks like and how to initialise the system to this state.

Each tutorial after this simulates an increasingly more complex atomic system with hyperfine structure.
