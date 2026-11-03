#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
#
# MIT License
#
# Copyright (c) 2026 Roger Cooperman and A.E. Petrunina V.V.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Lemma A.3 - adding an exotic test pair is an extra input."""
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

C = Checker("Lemma A.3 (Adding an exotic test pair is an extra input)",
            "Exercise 9 alone does not assert existence of an exotic partner for RP^4; Input B supplies it.",
            "Shows that every quantity Exercise 9 computes is homotopy invariant and identical for R and Q, "
            "so the exercise cannot produce or separate Q. " + INPUT_NOTE)
a = TruncPoly.gen(4)
w_R = (TruncPoly.one(4) + a) ** 5
# Q: computed from Wu's formula (homotopy invariant data), cf. Lemma 2.7
v = TruncPoly([1, 1, 1, 0, 0], 4)
w_Q = TruncPoly([0], 4)
for j, vj in enumerate(v.c):
    if vj:
        for i in range(j + 1):
            w_Q = w_Q + sq_power_of_a(j, i, 4)
C.check("total SW classes of R and Q agree", w_R == w_Q, repr(w_Q))
C.check("Exercise 9 outputs (orientable?, Pin+?, Pin-?) agree for R and Q",
        (w_R.deg(1), w_R.deg(2), (w_R.deg(2) + w_R.deg(1)) % 2) ==
        (w_Q.deg(1), w_Q.deg(2), (w_Q.deg(2) + w_Q.deg(1)) % 2))
C.check("the separating information (±1 vs ±9) comes only from Input B", R_CLASSES != Q_CLASSES)
C.done()
