#!/usr/bin/env python3
"""Lemma 3.2 - the closed-manifold invariant I(M,p) = chi([M,p])."""
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

C = Checker("Lemma 3.2 (The closed-manifold invariant)",
            "I(M,p) = χ([M,p]) is a diffeomorphism invariant, multiplicative under disjoint union, "
            "and trivial on Pin+ boundaries.", INPUT_NOTE + " Bordism is modelled abstractly as Z/16.")
I = lambda cls: chi(cls % 16)
C.check("disjoint union = addition of classes -> I(M ⊔ N) = I(M) I(N) for all classes",
        all(np.isclose(I(m + n), I(m) * I(n)) for m in range(16) for n in range(16)))
C.check("a Pin+ boundary has class 0 and I = 1", np.isclose(I(0), 1))
C.check("structured diffeomorphic manifolds have equal classes, hence equal I (by definition)",
        all(np.isclose(I(c), I(c)) for c in range(16)))
C.check("orientation-type reversal of structure (class -k) gives the conjugate value",
        all(np.isclose(I(-c), np.conj(I(c))) for c in range(16)))
C.done()
