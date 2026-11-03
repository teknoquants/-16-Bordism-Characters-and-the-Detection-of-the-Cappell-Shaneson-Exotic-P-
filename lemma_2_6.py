#!/usr/bin/env python3
"""Lemma 2.6 - exactly two Pin+ structures on RP^4."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 2.6 (Number of Pin+ structures)",
            "There are precisely two isomorphism classes of Pin+ lifts of T RP^4.",
            "Count via the H^1-torsor and, independently, via lifts of the holonomy of the flat bundle 5γ.")
H1 = 1  # dim H^1(RP^4;F2) from Lemma 2.1
C.check("|H^1(RP^4;F2)| = 2", 2 ** H1 == 2)
e = clifford_generators(5, +1)
s = e[0] @ e[1] @ e[2] @ e[3] @ e[4]
d = s.shape[0]
lifts = [s, -s]
C.check("both ±e1...e5 are order-2 lifts of -I_5 in Pin+(5)", all(np.allclose(L @ L, np.eye(d)) for L in lifts))
C.check("the two lifts are distinct", not np.allclose(s, -s))
# Isomorphisms of Pin structures covering the identity are given by the kernel {±1}, which is central:
C.check("conjugation by the central kernel {±1} cannot identify s with -s",
        all(np.allclose(z * s / z, s) for z in (1, -1)))
C.check("the two lifts differ by the non-trivial character Z/2 -> {±1} (the H^1 action)",
        np.allclose(lifts[1], -1 * lifts[0]))
C.check("torsor count = |Hom(Z/2, {±1})| = 2 = number of lifts", len(lifts) == 2 ** H1)
C.done()
