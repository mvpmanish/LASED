Using LASED from AI agents
==========================

LASED is designed to be driven by AI agents as well as people. The repository root contains
`AGENTS.md <https://github.com/mvpmanish/LASED/blob/main/AGENTS.md>`__, a one-page reference with
units, conventions, an API cheat sheet, recipes, a validation checklist and common pitfalls.

Key points:

* Units are ns, Grad/s (rad/ns), m, mW/mm², m/s and gauss.
* Results come back as plain NumPy arrays of shape ``(len(time), n, n)``.
* ``AtomicSystem.describe()`` returns a JSON-serialisable summary of the model, including every
  coupling.
* Always check Hermiticity, the trace and positivity of the result, and check convergence with
  respect to ``rtol``, ``n_ctrl`` and the ensemble size.
