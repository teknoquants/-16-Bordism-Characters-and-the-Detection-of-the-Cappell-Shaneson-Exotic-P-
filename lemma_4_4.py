#!/usr/bin/env python3
"""Lemma 4.4 - a scalar conjugation check is insufficient."""
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

C = Checker("Lemma 4.4 (A scalar conjugation check is insufficient)",
            "For k in Z/16, conj(χ(k)) = χ(-k); this alone does not define T on state spaces.")
C.check("conj(χ(k)) = χ(-k) for all k", all(np.isclose(np.conj(chi(k)), chi(-k % 16)) for k in range(16)))
C.check("conjugation sends the theory χ to the different theory χ^-1 = χ^15",
        not np.allclose([chi(k) for k in range(16)], [np.conj(chi(k)) for k in range(16)]))
C.note("the identity involves only closed 4-manifold values (no 3-manifold state-space data)", "structural observation")
C.done()
