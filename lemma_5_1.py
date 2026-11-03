#!/usr/bin/env python3
"""Lemma 5.1 - J1^2 = -I and J1 M J1^-1 = Y conj(M) Y."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.1 (Square and conjugation formula)",
            "J1^2 = -I, and J1 M J1^{-1} = Y conj(M) Y for every complex 2x2 matrix M.")
J1 = Antilinear(Y)
C.check("conj(Y) = -Y and Y^2 = I", np.allclose(np.conj(Y), -Y) and np.allclose(Y @ Y, I2))
C.check("J1^2 = Y conj(Y) = -I", np.allclose(J1.square(), -I2))
J1inv = Antilinear(-Y)                       # (YK)^{-1} = KY = conj(Y) K = -Y K
C.check("(YK)^{-1} = KY = -YK", np.allclose(J1 @ J1inv, I2) and np.allclose(J1inv @ J1, I2))
ok = True
for _ in range(500):
    M = random_complex((2, 2))
    conjugated = J1 @ Antilinear(M @ (-Y))   # J1 ∘ M ∘ J1^{-1}  (M ∘ (-Y K) = (-M Y) K)
    ok &= np.allclose(conjugated, Y @ np.conj(M) @ Y)
C.check("J1 M J1^{-1} = Y conj(M) Y for 500 random M", ok)
v = random_complex(2)
C.check("antilinearity: J1(i v) = -i J1(v)", np.allclose(J1(1j * v), -1j * J1(v)))
C.done()
