#!/usr/bin/env python3
"""Lemma 3.1 - chi(k) = zeta^k is a well-defined injective character of Z/16."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 3.1 (An explicit finite character)",
            "With ζ = e^{2πi/16}, χ(k) = ζ^k is a well-defined injective homomorphism Z/16 -> C^×.")
z = np.exp(2j * np.pi / 16)
chi = lambda k: z ** k
C.check("well-defined: χ(k + 16m) = χ(k) for k in 0..15, m in -3..3",
        all(np.isclose(chi(k + 16 * m), chi(k)) for k in range(16) for m in range(-3, 4)))
C.check("multiplicative: χ(k + l) = χ(k)χ(l) for all k, l",
        all(np.isclose(chi((k + l) % 16), chi(k) * chi(l)) for k in range(16) for l in range(16)))
C.check("kernel is {0}: χ(k) = 1 iff 16 | k", [k for k in range(16) if np.isclose(chi(k), 1)] == [0])
C.check("ζ is a primitive 16th root of unity", min(k for k in range(1, 17) if np.isclose(z ** k, 1)) == 16)
C.done()
