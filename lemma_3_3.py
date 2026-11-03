#!/usr/bin/env python3
"""Lemma 3.3 - separation independent of labelling and of the identification."""
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

C = Checker("Lemma 3.3 (Separation independent of labelling)",
            "No choice of Pin+ structure on R gives the same value of I as any choice on Q.", INPUT_NOTE)
C.check("exponent sets are {1,15} and {7,9}", R_CLASSES == {1, 15} and Q_CLASSES == {7, 9})
C.check("they are disjoint", R_CLASSES.isdisjoint(Q_CLASSES))
vals = lambda S, u=1: {np.round(chi(k, u), 12) for k in S}
C.check("value sets {ζ, ζ^-1} and {ζ^9, ζ^-9} are disjoint", vals(R_CLASSES).isdisjoint(vals(Q_CLASSES)))
units = [u for u in range(16) if np.gcd(u, 16) == 1]
C.check("disjoint under every automorphism of Z/16 (u odd)",
        all({u * k % 16 for k in R_CLASSES}.isdisjoint({u * k % 16 for k in Q_CLASSES}) for u in units))
C.check("disjoint for every injective character χ^u (u odd)",
        all(vals(R_CLASSES, u).isdisjoint(vals(Q_CLASSES, u)) for u in units))
bad = [u for u in range(16) if not vals(R_CLASSES, u).isdisjoint(vals(Q_CLASSES, u))]
C.check("control: exactly the even u (characters killing 8) fail to separate", bad == list(range(0, 16, 2)),
        f"failing u = {bad}")
C.check("all 4 pairings (p on R, q on Q) give different values",
        all(not np.isclose(chi(p), chi(q)) for p in R_CLASSES for q in Q_CLASSES))
C.done()
