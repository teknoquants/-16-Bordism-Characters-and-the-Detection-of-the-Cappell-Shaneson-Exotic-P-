#!/usr/bin/env python3
"""Lemma A.4 - from a test object to an existence question."""
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

C = Checker("Lemma A.4 (From a test object to an existence question)",
            "The exercise, the symmetry sign and an exotic pair motivate, but do not prove, the existence assertion.",
            "Illustrative: the three ingredients are all compatible with theories that separate R and Q "
            "and with theories that do not, so they cannot by themselves force existence. " + INPUT_NOTE)
separates = lambda u: all(not np.isclose(chi(p, u), chi(q, u)) for p in R_CLASSES for q in Q_CLASSES)
C.check("a separating Pin+ character exists (u = 1)", separates(1))
C.check("non-separating Pin+ characters also exist (u = 0, 2)", not separates(0) and not separates(2))
C.note("so the ingredients identify what to search for; a witness (Input E + Lemma 3.3) is extra")
C.done()
