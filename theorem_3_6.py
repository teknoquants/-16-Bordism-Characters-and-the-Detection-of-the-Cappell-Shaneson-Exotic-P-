#!/usr/bin/env python3
"""Theorem 3.6 - Pin+ detection as a consequence of established results."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

# ---- External Input B (cited, not proved here): Stolz [6] / Debray [7, Ex. 3.4] / Kirby-Taylor [8]
OMEGA = 16                      # Ω_4^{Pin+} ≅ Z/16
R_CLASSES = {1, -1 % 16}        # the two Pin+ structures on RP^4
Q_CLASSES = {9, -9 % 16}        # the two Pin+ structures on the Cappell-Shaneson RP^4
INPUT_NOTE = ("Input B (Ω4^Pin+ ≅ Z/16; RP^4 ↦ ±1, Cappell-Shaneson RP^4 ↦ ±9; homeomorphism) "
              "is taken from the literature and is NOT validated by this script.")
zeta = roots16(1)
def chi(k, u=1):
    return roots16(u * k)

C = Checker("Theorem 3.6 (Pin+ detection as a consequence of established results)",
            "There is an invertible 4D Pin+ super TQFT whose closed values distinguish RP^4 from the "
            "Cappell-Shaneson RP^4 for every pairing of Pin+ structures; its closed values are 16th roots of unity.",
            INPUT_NOTE + " Input E (existence of the functor) is external; the categorical/time-reversal "
            "clauses of Problem 5 are explicitly NOT claimed.")
pairs = [(p, q) for p in sorted(R_CLASSES) for q in sorted(Q_CLASSES)]
C.check("Lemma 2.6/2.7: each manifold has exactly two Pin+ structures", len(R_CLASSES) == len(Q_CLASSES) == 2)
C.check("Lemma 3.3: all 4 pairings separated", all(not np.isclose(chi(p), chi(q)) for p, q in pairs),
        str(pairs))
C.check("Lemma 3.1: closed values lie in μ16",
        all(np.isclose(chi(k) ** 16, 1) for k in R_CLASSES | Q_CLASSES))
C.check("the separating values are not all equal to 1 (the theory is non-trivial)",
        any(not np.isclose(chi(k), 1) for k in R_CLASSES | Q_CLASSES))
C.note("homeomorphism of R and Q", "Input B(i) - Hambleton-Kreck-Teichner; external, not validated")
C.done()
