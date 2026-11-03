#!/usr/bin/env python3
"""Lemma 1.1 - T^2 = P is not an involution when the odd sector is nonzero."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 1.1 (The wording 'involution')",
            "On a space with a nonzero odd sector, a map satisfying T^2 = P is not an "
            "involution in the literal sense T^2 = I.")

# Super vector space C^2 (even) + C^2 (odd); P = diag(1,1,-1,-1).
P = np.diag([1, 1, -1, -1]).astype(complex)
U = np.zeros((4, 4), complex)
U[:2, :2] = I2          # even sector: plain conjugation K, squares to +1
U[2:, 2:] = Y           # odd sector: Y K, squares to -1
T = Antilinear(U)
T2 = T.square()
C.check("T is antilinear with T^2 = P", np.allclose(T2, P))
C.check("T^2 != I", not np.allclose(T2, np.eye(4)))
v = np.array([0, 0, 1, 0], complex)
C.check("for an odd vector v, T^2 v = -v != v", np.allclose(T2 @ v, -v) and not np.allclose(T2 @ v, v))
# If T^2 = I and T^2 = P simultaneously, P = I, i.e. no odd sector.
for odd_dim in range(0, 4):
    Pk = np.diag([1] * (4 - odd_dim) + [-1] * odd_dim)
    C.check(f"P = I  <=>  odd dimension 0  (odd dim {odd_dim})",
            np.allclose(Pk, np.eye(4)) == (odd_dim == 0))
C.done()
