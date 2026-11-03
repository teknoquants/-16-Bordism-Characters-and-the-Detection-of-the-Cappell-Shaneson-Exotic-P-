#!/usr/bin/env python3
"""Lemma 4.2 - the Pin+ obstruction data does not construct or single out a TQFT."""
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

C = Checker("Lemma 4.2 (What the Pin+ obstruction does not show)",
            "w2 = 0 and w2 + w1^2 != 0 do not construct a TQFT or certify its finite categorical origin.",
            "Illustrative: the same obstruction data is compatible with 16 different invertible theories "
            "χ^u, which differ on RP^4 and include the trivial one. " + INPUT_NOTE)
theories = {u: tuple(np.round(chi(k, u), 10) for k in range(16)) for u in range(16)}
C.check("16 pairwise different character theories share the same Pin+ admissibility data",
        len(set(theories.values())) == 16)
C.check("the trivial theory (u = 0) is among them and does not separate R from Q",
        all(np.isclose(chi(k, 0), 1) for k in R_CLASSES | Q_CLASSES))
C.check("so (w1, w2) alone cannot determine the theory's values on RP^4",
        len({np.round(chi(1, u), 10) for u in range(16)}) == 16)
C.done()
