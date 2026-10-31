Physics and conventions
=======================

Equations of motion
-------------------

LASED solves the Lindblad master equation for the density matrix :math:`\rho` of all sub-states,

.. math::

   \dot\rho = -i[H(t), \rho] + \sum_k \left( C_k \rho C_k^\dagger - \tfrac12\{C_k^\dagger C_k, \rho\} \right)
              - \tfrac12 \{K, \rho\},

in units where frequencies are in Grad/s (rad/ns) and times in ns.

* :math:`H` is the Hamiltonian in the rotating frame of the lasers, made Hermitian by construction:

  .. math::

     H = \sum_i (\omega_i - f_i)|i\rangle\langle i| + \sum_{\rm lasers} \left( V + V^\dagger \right),
     \qquad V_{ul} = \frac{\Omega_0}{2} A(t) \sum_q a_q\, c(u, l, q).

  Here :math:`f_i` is the frame frequency of state :math:`i` (sums of laser frequencies along the
  laser graph), :math:`A(t)` is the complex pulse envelope, :math:`a_q` are the polarisation
  amplitudes and :math:`c(u,l,q)` = ``coupling(u, l, q)`` is the dipole matrix element from the
  Wigner–Eckart theorem (3j and 6j symbols).
* :math:`C_k` are spontaneous-emission jump operators,
  :math:`C_q = \sqrt{\Gamma}\sum_{e,g} c(e,g,q)/\sqrt{N_e}\, |g\rangle\langle e|`, with
  :math:`N_e = \sum_{g,q} c(e,g,q)^2`, plus dephasing operators.
* :math:`K` collects losses out of the modelled system (``tau_f``, ``tau_b``, ``Loss``).

The superoperator is assembled in row-major vectorisation, :math:`\mathrm{vec}(\rho)_{(i-1)n+(j-1)} = \rho_{ij}`
(the same convention as ``index(i, j, n)``), as
:math:`-i(H\otimes 1 - 1\otimes H^T) + \sum_k [C_k\otimes C_k^* - \tfrac12 (C_k^\dagger C_k\otimes 1 + 1\otimes (C_k^\dagger C_k)^T)]`.

Polarisation
------------

The quantisation axis is :math:`z`. The polarisation vector is expanded as
:math:`\boldsymbol\epsilon = \sum_q a_q \mathbf u_q`, with :math:`\mathbf u_{\pm1} = (\hat x \pm i\hat y)/\sqrt2`
and :math:`\mathbf u_0 = \hat z`, so :math:`a_q = \mathbf u_q^*\cdot\boldsymbol\epsilon`.
:math:`a_q` drives :math:`\Delta m = m_u - m_l = q`. Use :func:`LASED.polarisation.jones`,
:func:`~LASED.polarisation.linear`, :func:`~LASED.polarisation.circular` and
:func:`~LASED.polarisation.elliptical` to build them. The legacy arguments ``Q`` and ``rabi_factors``
are the same amplitudes, ``{Q[i]: rabi_factors[i]}``.

This convention is verified in the test suite in two ways:

* against an independent QuTiP model built from Cartesian dipole operators, for elliptical and
  general complex polarisations;
* through rotation covariance: rotating the density matrix with Wigner D-matrices and driving with π
  light reproduces direct driving with the rotated linear polarisation.

Elliptical light: what was wrong in v1
--------------------------------------

LASED v1 wrote the equations out term by term and used the coupling :math:`\Omega_{eg}` where the
complex conjugate :math:`\Omega_{eg}^*` belongs (the :math:`H_{ge}` elements). For real amplitudes
the two are equal. For complex amplitudes, as with elliptical light or rotated linear polarisation,
the generator stops preserving Hermiticity. The trace stays equal to 1, but individual populations
grow without bound and become negative. Version 2 builds :math:`H` as a Hermitian matrix and derives
the equations from it, so this cannot happen.

A second, smaller bug was also fixed. Spontaneous emission did not transfer coherences between ground
sub-states whose excited partners differ in :math:`m`, which happens whenever the light mixes
:math:`\sigma^+` and :math:`\sigma^-`. Lindblad jump operators include this automatically.

Why one decay operator per photon polarisation
----------------------------------------------

A tempting shortcut is one jump operator per decay path, :math:`L_{eg} = \sqrt{\Gamma_{eg}}\,|g\rangle\langle e|`.
It gives the right populations, but it is wrong whenever the levels have Zeeman (or hyperfine)
structure. A photon emitted with polarisation :math:`q` does not record which of several degenerate
paths produced it, so those paths interfere. The correct operators are one per polarisation,

.. math::

   C_q = \sqrt{\Gamma} \sum_{e,g} \frac{c(e,g,q)}{\sqrt{N_e}} |g\rangle\langle e|, \qquad q \in \{-1, 0, +1\}.

These transfer coherences between excited sub-states to the ground sub-states. Example: an excited
superposition :math:`(|e,-\tfrac12\rangle + |e,+\tfrac12\rangle)/\sqrt2` of a :math:`J=\tfrac12\to\tfrac12`
transition leaves a ground-state coherence :math:`|\rho_{g-,g+}| = 1/6` after it decays. Per-path
operators predict zero. ``tests/test_backend_storage.py`` checks both. LASED v1 had a version of this
error for light that mixes :math:`\sigma^+` and :math:`\sigma^-`.

How the operators are stored and built
--------------------------------------

The user-facing objects are unchanged: ``State`` objects declared with keyword quantum numbers,
``generateSubStates``, ``LaserAtomSystem``. Internally, an ``AtomicSystem`` converts the states once into
a :class:`LASED.basis.Basis`, with one integer array per (doubled) quantum number. Dipole coupling blocks
for every pair of manifolds are then computed in one vectorised call to the Racah formulas
(:mod:`LASED.wigner`), and :math:`H`, :math:`C_q` and :math:`K` are filled by array indexing. The
Liouvillian is assembled from them with sparse Kronecker products.

Rotating frame, detuning and Doppler shifts
-------------------------------------------

The frame frequencies come from a breadth-first search over the graph of lasers. Every laser edge
lower→upper adds that laser's frequency. Positive ``detuning`` raises the laser frequency. An atom
with velocity :math:`\mathbf v` sees each laser shifted by :math:`-\mathbf k\cdot\mathbf v`. A closed
loop of lasers has no time-independent frame and is rejected.

Further reading
---------------

* D. A. Steck, *Quantum and Atom Optics* (lecture notes), and *Rubidium/Caesium D Line Data*: spherical
  basis, dipole matrix elements and the optical Bloch equations.
* M. Auzinsh, D. Budker and S. M. Rochester, *Optically Polarized Atoms* (Oxford, 2010): polarisation,
  Zeeman coherences and multi-level density matrices.
* H.-P. Breuer and F. Petruccione, *The Theory of Open Quantum Systems* (Oxford, 2002): the Lindblad form
  and why it preserves trace and positivity.
* D. Manzano, "A short introduction to the Lindblad master equation", AIP Advances **10**, 025106 (2020).
* R. Loudon, *The Quantum Theory of Light* (Oxford): the interaction Hamiltonian and spontaneous emission.
* J. R. Johansson, P. D. Nation and F. Nori, "QuTiP 2: A Python framework for the dynamics of open quantum
  systems", Comput. Phys. Commun. **184**, 1234 (2013).
* S. Stenholm, "The semiclassical theory of laser cooling", Rev. Mod. Phys. **58**, 699 (1986).
