#!/usr/bin/env python3
"""Lemma 3.4 - averaged detector A(R) = cos(pi/8), A(Q) = -cos(pi/8)."""
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

C = Checker("Lemma 3.4 (A real unstructured numerical separation)",
            "A(R) = cos(π/8), A(Q) = -cos(π/8); for any identification/injective character the values "
            "are ±cos(πu/8), u odd, and are unequal.", INPUT_NOTE)
A = lambda S, u=1: 0.5 * sum(chi(k, u) for k in S)
C.check("A(R) = cos(π/8) (real)", np.isclose(A(R_CLASSES), np.cos(np.pi / 8)))
C.check("A(Q) = -cos(π/8) (real)", np.isclose(A(Q_CLASSES), -np.cos(np.pi / 8)))
C.check("ζ^9 = -ζ and ζ^-9 = -ζ^-1", np.isclose(zeta ** 9, -zeta) and np.isclose(zeta ** -9, -1 / zeta))
C.check("cos(π/8) > 0", np.cos(np.pi / 8) > 0)
for u in range(1, 16, 2):
    aR, aQ = A(R_CLASSES, u), A(Q_CLASSES, u)
    C.check(f"u={u}: A(R)=cos(πu/8), A(Q)=-A(R), and A(R) != 0",
            np.isclose(aR, np.cos(np.pi * u / 8)) and np.isclose(aQ, -aR) and abs(aR) > 1e-6)
C.check("independent of labelling: permuting the two lifts leaves the sums fixed",
        np.isclose(A(sorted(R_CLASSES)), A(sorted(R_CLASSES, reverse=True))))
C.done()
