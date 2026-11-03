#!/usr/bin/env python3
"""Lemma B.1 - internal arithmetic checks."""
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

C = Checker("Lemma B.1 (Internal arithmetic checks)",
            "The characteristic-class expansion, the exponent comparison and the averaged detector are consistent.")
C.check("binomial coefficients of (1+a)^5 mod 2 are (1,1,0,0,1,1)", [comb(5, k) % 2 for k in range(6)] == [1, 1, 0, 0, 1, 1])
C.check("-1 ≡ 15 and -9 ≡ 7 (mod 16)", (-1) % 16 == 15 and (-9) % 16 == 7)
C.check("neither 7 nor 9 equals 1 or 15", {7, 9}.isdisjoint({1, 15}))
C.check("adding 8 to an exponent multiplies the root of unity by -1",
        all(np.isclose(roots16(k + 8), -roots16(k)) for k in range(16)))
C.check("hence A(Q) = -A(R)", np.isclose(0.5 * (roots16(9) + roots16(7)), -0.5 * (roots16(1) + roots16(15))))
C.done()
