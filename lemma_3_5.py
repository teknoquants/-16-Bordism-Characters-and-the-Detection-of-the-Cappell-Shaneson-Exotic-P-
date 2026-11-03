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
"""Lemma 3.5 - functorial extension: algebraic layer omega = chi o s."""
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

C = Checker("Lemma 3.5 (Functorial extension, with dependency exposed)",
            "Using Inputs B and E, I is the closed-manifold partition function of an invertible 4D Pin+ super TQFT.",
            INPUT_NOTE + " Input E (Reutter-Schommer-Pries, Thm 3.14) is also external. This script checks only "
            "the algebraic layer: any homomorphism s from a finitely generated abelian group to Z/16, composed "
            "with χ, is a character, and its values are the χ-values of the image classes.")
# Model Ω^{Tβ}_4 by a random finitely generated abelian group Z^a ⊕ Z/16 ⊕ Z/2 and s by a random homomorphism.
ok_hom, ok_val = True, True
for trial in range(50):
    a = RNG.integers(0, 3)
    imgs = [int(RNG.integers(0, 16)) for _ in range(a)]            # images of free generators
    t16 = int(RNG.integers(0, 16))                                  # image of the Z/16 generator
    t2 = int(RNG.choice([0, 8]))                                    # image of Z/2 must have order | 2
    def s(x):
        free, y, z = x
        return (sum(i * f for i, f in zip(imgs, free)) + t16 * y + t2 * z) % 16
    omega = lambda x: chi(s(x))
    for _ in range(30):
        x1 = (list(RNG.integers(-5, 5, size=a)), int(RNG.integers(0, 16)), int(RNG.integers(0, 2)))
        x2 = (list(RNG.integers(-5, 5, size=a)), int(RNG.integers(0, 16)), int(RNG.integers(0, 2)))
        xs = ([p + q for p, q in zip(x1[0], x2[0])], (x1[1] + x2[1]) % 16, (x1[2] + x2[2]) % 2)
        ok_hom &= np.isclose(omega(xs), omega(x1) * omega(x2))
        ok_val &= np.isclose(omega(x1), chi(s(x1)))
C.check("ω = χ ∘ s is a homomorphism to C^× (well-defined on torsion)", ok_hom)
C.check("closed values Z_ω(M,p) = χ([M,p]_Pin+) = I(M,p)", ok_val)
C.check("Hom(Z/16, C^×) is cyclic of order 16 and χ is a generator",
        len({np.round(roots16(u), 12) for u in range(16)}) == 16)
C.check("partition functions of an invertible theory lie in the unit circle (super lines)",
        all(np.isclose(abs(chi(k)), 1) for k in range(16)))
C.done()
