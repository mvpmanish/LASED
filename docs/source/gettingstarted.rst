================
Getting Started
================

Getting Started
================

What is LASED?
---------------

LASED stands for Laser-Atom interaction Simulator using quantum ElectroDynamics. LASED is a python library which can:

* Calculate the time evolution an atomic system interacting with a laser.
* Generate the equations of motion of an atom-laser system.
* Rotate an atomic system to a different reference frame.

LASED can simulate any atomic system. The sub-states, angular momenta, spin, and energies are provided by the user to simulate the atomic system.

Laser parameters also need to be specified to build the laser-atom system such as detuning from the transition frequency, polarisation, and laser intensity/power.

LASED can simulate a Gaussian beam profile and doppler averaging of the atoms.

Installation
--------------

You can easily install LASED by opening up a terminal and running::

  pip install LASED

The source code can be viewed `here <https://github.com/mvpmanish/LASED>`__

Using LASED
------------

Start by going to :doc:`tutorials`.

The first tutorial is a guide on how to simulate one of the simplest excitations of an atom with a laser and shows how to use LASED.

The second tutorial shows how average over the laser beam and doppler profile of the atoms. It also introduces rotating the system to different reference frames using the Wigner rotation matrix.

The third tutorial simulates a metastable Helium atom with some real data on what the initial state looks like and how to initialise the system to this state.

Each tutorial after this simulates an increasingly more complex atomic system with hyperfine structure.
