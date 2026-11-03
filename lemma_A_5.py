#!/usr/bin/env python3
"""Lemma A.5 - a finite list of agreeing invariants cannot close the gap."""
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

C = Checker("Lemma A.5 (The finite numerical genesis script cannot close the gap)",
            "Agreement of a chosen finite list of invariants on an exotic pair does not prove agreement "
            "for every TQFT in Problem 5.",
            "Builds a finite list of homotopy-invariant quantities (Euler characteristic, all SW numbers, "
            "Pin± admissibility, untwisted Dijkgraaf-Witten counts |Hom(π1,G)|/|G| for four finite groups) "
            "that agree on R and Q, then exhibits the Pin+ bordism character that does not. " + INPUT_NOTE)
a = TruncPoly.gen(4)
w = (TruncPoly.one(4) + a) ** 5
def sw_numbers(w):
    w1, w2, w3, w4 = (w.deg(k) for k in (1, 2, 3, 4))
    # evaluate degree-4 monomials on the fundamental class (H^4 = F2, a^4 -> 1)
    return {"w1^4": w1, "w1^2 w2": (w1 * w2) % 2, "w2^2": w2, "w1 w3": (w1 * w3) % 2, "w4": w4}
inv_R = {"euler": 1, **sw_numbers(w)}
inv_Q = {"euler": 1, **sw_numbers(w)}         # homotopy invariants (Lemma 2.7)
def count_hom_Z2(G_elems, mul, e):
    return sum(1 for g in G_elems if mul(g, g) == e)      # Hom(Z/2, G) = elements with g^2 = 1
groups = {
    "Z/2": (list(range(2)), lambda x, y: (x + y) % 2, 0),
    "Z/3": (list(range(3)), lambda x, y: (x + y) % 3, 0),
    "Z/4": (list(range(4)), lambda x, y: (x + y) % 4, 0),
    "S3": (list(itertools.permutations(range(3))), lambda p, q: tuple(p[q[i]] for i in range(3)), (0, 1, 2)),
}
for name, (G, mul, e) in groups.items():
    dw = count_hom_Z2(G, mul, e) / len(G)
    inv_R[f"DW[{name}]"] = dw
    inv_Q[f"DW[{name}]"] = dw                  # π1 = Z/2 for both
inv_R["Pin+ admissible"] = inv_Q["Pin+ admissible"] = (w.deg(2) == 0)
C.check(f"all {len(inv_R)} listed invariants agree on R and Q", inv_R == inv_Q,
        ", ".join(f"{k}={v}" for k, v in inv_R.items()))
C.check("yet the Pin+ bordism character separates them (Lemma 3.3)",
        all(not np.isclose(chi(p), chi(q)) for p in R_CLASSES for q in Q_CLASSES))
C.done()
