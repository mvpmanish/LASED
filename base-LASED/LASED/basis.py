'''
Internal storage of the atomic basis as arrays.

Users declare sub-states with keyword-argument ``State`` objects (like reading them off a Grotrian
diagram). Internally LASED converts the list once into a ``Basis``: one NumPy array per quantum
number, so that Hamiltonians, dipole matrices and collapse operators are built with array operations
instead of Python loops over matrix elements.

Quantum numbers are stored doubled (``L2 = 2L`` etc.) as integers so half-integers are exact.
'''

from dataclasses import dataclass
import numpy as np

from LASED.wigner import dipoleCoupling

QUANTUM_NUMBERS = ("L", "S", "J", "I", "F", "m")


def _twice(x):
    r = np.rint(2*np.asarray(x, dtype=float)).astype(np.int64)
    if np.any(np.abs(r - 2*np.asarray(x, dtype=float)) > 1e-9):
        raise ValueError(f"Quantum numbers must be integers or half-integers, got {x}")
    return r


@dataclass
class Basis:
    """Structure-of-arrays view of a list of ``State`` objects, ordered by label.

    Attributes:
        states (list of State): The user's State objects, sorted by label.
        label (ndarray[int]): Labels (1..n).
        w (ndarray[float]): Angular frequencies in Grad/s.
        qn2 (dict[str, ndarray[int]]): Doubled quantum numbers "L", "S", "J", "I", "F", "m".
    """
    states: list
    label: np.ndarray
    w: np.ndarray
    qn2: dict

    @classmethod
    def fromStates(cls, states):
        states = sorted(states, key=lambda s: s.label)
        label = np.array([s.label for s in states], dtype=np.int64)
        if not np.array_equal(label, np.arange(1, len(states) + 1)):
            raise ValueError(f"State labels must be exactly 1..n with no gaps or duplicates, got {label.tolist()}.")
        w = np.array([s.w for s in states], dtype=float)
        qn2 = {k: _twice([getattr(s, k) for s in states]) for k in QUANTUM_NUMBERS}
        return cls(states, label, w, qn2)

    @property
    def n(self):
        return len(self.label)

    def indices(self, states):
        """Matrix indices (label - 1) of a list of State objects (or labels)."""
        return np.array([getattr(s, "label", s) - 1 for s in states], dtype=np.int64)

    def subset(self, idx):
        """Doubled quantum numbers of the sub-states at matrix indices ``idx``."""
        return {k: v[idx] for k, v in self.qn2.items()}

    def dipoleBlock(self, upper, lower, q):
        """Coupling coefficients coupling(u, l, q) for all u in ``upper`` and l in ``lower``.

        Returns:
            ndarray: shape (len(upper), len(lower)).
        """
        return dipoleCoupling(self.subset(self.indices(upper)), self.subset(self.indices(lower)), q)

    def table(self):
        """The basis as a table: a pandas DataFrame if pandas is installed, else a dict of arrays."""
        data = {"label": self.label, "w": self.w}
        data.update({k: self.qn2[k]/2 for k in QUANTUM_NUMBERS})
        try:
            import pandas as pd
            return pd.DataFrame(data).set_index("label")
        except ImportError:
            return data
