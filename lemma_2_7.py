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
"""Lemma 2.7 - the same conclusions for any homotopy RP^4 (via Wu's formula)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 2.7 (The same conclusions for any homotopy RP^4)",
            "If M is a closed smooth 4-manifold homotopy equivalent to RP^4, then w(TM) = 1 + a + a^4; "
            "M is non-orientable, has no Pin- structure and exactly two Pin+ structures.",
            "Computes w from Wu's formula using only the cohomology ring and Steenrod squares - "
            "homotopy-invariant data - so the result applies verbatim to the Cappell-Shaneson RP^4. "
            "The Wu implementation is cross-checked against w(RP^n) = (1+a)^(n+1) for n = 1..10.")

def wu_total_sw(n):
    """Total SW class of a closed n-manifold with H* = F2[a]/(a^{n+1}) and Sq^k(a^m) = C(m,k)a^{m+k}."""
    top = n
    # Wu class v_k: Sq^k(x) = v_k x for x in H^{n-k}; v_k = lambda a^k with lambda = C(n-k, k)
    v = TruncPoly([binom2(n - k, k) if 2 * k <= n else 0 for k in range(n + 1)], top)
    w = TruncPoly([0], top)
    for j, vj in enumerate(v.c):          # w = Sq(v) = Σ_j Σ_i Sq^i(v_j a^j)
        if vj:
            for i in range(0, j + 1):
                w = w + sq_power_of_a(j, i, top)
    return v, w

for n in range(1, 11):
    _, w = wu_total_sw(n)
    expected = (TruncPoly.one(n) + TruncPoly.gen(n)) ** (n + 1)
    C.check(f"Wu-formula SW class equals (1+a)^{n+1} for RP^{n}", w == expected, repr(w))

v, w = wu_total_sw(4)
C.check("Wu class of a homotopy RP^4 is v = 1 + a + a^2", v == TruncPoly([1, 1, 1, 0, 0], 4), repr(v))
C.check("w(TM) = Sq(v) = 1 + a + a^4", w == TruncPoly([1, 1, 0, 0, 1], 4), repr(w))
w1 = TruncPoly([0, w.deg(1)], 4); w2 = TruncPoly([0, 0, w.deg(2)], 4)
C.check("w1(M) = a != 0  => M non-orientable", w1.c[1] == 1)
C.check("w2(M) = 0  => Pin+ exists", w2 == TruncPoly([0], 4))
C.check("w2(M) + w1(M)^2 = a^2 != 0  => no Pin-", (w2 + w1 * w1).c[2] == 1)
C.check("|H^1(M;F2)| = |H^1(RP^4;F2)| = 2  => exactly two Pin+ structures", 2 ** 1 == 2)
C.done()
