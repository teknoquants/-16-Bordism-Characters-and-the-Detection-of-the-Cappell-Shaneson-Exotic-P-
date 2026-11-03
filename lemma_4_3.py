#!/usr/bin/env python3
"""Lemma 4.3 - three meanings of finiteness."""
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

C = Checker("Lemma 4.3 (Three meanings of finiteness)",
            "Finite image of closed partition functions and finite-dimensional state spaces do not, by "
            "definition alone, verify a finite local categorical presentation.",
            "Illustrative: every character χ^u has finite image and 1-dimensional state spaces, so these "
            "two properties cannot distinguish theories with and without a presentation; a presentation "
            "is extra data that no scalar computation supplies. Not a theorem-proving script.")
for u in (1, 2, 3, 8):
    img = {np.round(chi(k, u), 10) for k in range(16)}
    C.check(f"χ^{u}: finite image of size {len(img)} (divides 16)", 16 % len(img) == 0)
C.note("invertible theories have 1-dimensional (super-line) state spaces in every case", "structural property of invertible theories, not computed")
C.check("the two properties hold uniformly in u, so they carry no information about a local presentation",
        all(16 % len({np.round(chi(k, u), 10) for k in range(16)}) == 0 for u in range(16)))
C.done()
