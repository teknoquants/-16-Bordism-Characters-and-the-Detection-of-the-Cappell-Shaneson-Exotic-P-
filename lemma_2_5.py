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
"""Lemma 2.5 - RP^4 is non-orientable, Pin+ but not Pin-."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 2.5 (All three conclusions of Exercise 9)",
            "RP^4 is non-orientable, admits Pin+ structures, and admits no Pin- structure.",
            "Two independent routes: (A) characteristic classes; (B) the flat bundle 5γ = TRP^4 ⊕ ε "
            "has holonomy -I_5 around the generator of π1 = Z/2, so a Pin± structure on it is a lift "
            "of -I_5 to an element of order 2 in Pin±(5). Pin structures are stable under adding ε.")
a = TruncPoly.gen(4)
w = (TruncPoly.one(4) + a) ** 5
w1, w2 = TruncPoly([0, w.deg(1)], 4), TruncPoly([0, 0, w.deg(2)], 4)
C.check("(A) w1 = a != 0  => non-orientable", w1.c[1] == 1)
C.check("(A) w2 = 0  => Pin+ exists", all(x == 0 for x in w2.c))
C.check("(A) w2 + w1^2 = a^2 != 0  => no Pin-", (w2 + w1 * w1) == TruncPoly([0, 0, 1], 4))

# (B) holonomy representation of Z/2 on R^5 is -I
C.check("(B) det(-I_5) = -1 => orientation character non-trivial", np.isclose(np.linalg.det(-np.eye(5)), -1))
for sign, name in ((+1, "Pin+"), (-1, "Pin-")):
    e = clifford_generators(5, sign)
    s = e[0] @ e[1] @ e[2] @ e[3] @ e[4]          # lifts the product of 5 reflections = -I_5
    dim = s.shape[0]
    sq = s @ s
    if sign == +1:
        C.check(f"(B) {name}: lift s = e1...e5 satisfies s^2 = +1 -> homomorphism Z/2 -> Pin+(5) exists",
                np.allclose(sq, np.eye(dim)))
    else:
        C.check(f"(B) {name}: both lifts ±e1...e5 square to -1 -> no lift of the holonomy",
                np.allclose(sq, -np.eye(dim)) and np.allclose((-s) @ (-s), -np.eye(dim)))
    # cocycle restricted to the diagonal element δ = (1,1,1,1,1): c(δ,δ) = C(5,2) + eps*5
    eps = 0 if sign == 1 else 1
    c_dd = (comb(5, 2) + eps * 5) % 2
    C.check(f"(B) {name}: pulled-back extension class on Z/2 is c(δ,δ) = {c_dd}",
            c_dd == eps)
    # On Z/2 the 2-cocycle with c(g,g)=1 is not a coboundary
    cz = lambda g, h, v=c_dd: v if (g == 1 and h == 1) else 0
    C.check(f"(B) {name}: that class is {'non-zero' if eps else 'zero'} in H^2(Z/2;F2)",
            is_2coboundary(cz, 1) == (eps == 0))
C.note("(B) H^2(BZ/2;F2) -> H^2(RP^4;F2) is an isomorphism (RP^4 is the 4-skeleton of RP^∞)", "cellular: both have one 2-cell and zero mod-2 coboundaries (see Lemma 2.1)")
C.done()
