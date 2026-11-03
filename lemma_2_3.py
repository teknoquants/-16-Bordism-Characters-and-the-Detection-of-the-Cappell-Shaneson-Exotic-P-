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
"""Lemma 2.3 - w(RP^4) = 1 + a + a^4."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 2.3 (Characteristic-class calculation)",
            "w(RP^4) = 1 + a + a^4; w1 = a, w2 = 0, w3 = 0, w4 = a^4.")
coeffs = [comb(5, k) for k in range(6)]
C.check("binomial coefficients of (1+a)^5 are 1,5,10,10,5,1", coeffs == [1, 5, 10, 10, 5, 1])
C.check("mod 2 they are 1,1,0,0,1,1", [c % 2 for c in coeffs] == [1, 1, 0, 0, 1, 1])
a = TruncPoly.gen(4)
w = (TruncPoly.one(4) + a) ** 5
C.check("(1+a)^5 = 1 + a + a^4 in F2[a]/(a^5)", w == TruncPoly([1, 1, 0, 0, 1], 4), repr(w))
C.check("w1 = a", w.deg(1) == 1)
C.check("w2 = 0", w.deg(2) == 0)
C.check("w3 = 0", w.deg(3) == 0)
C.check("w4 = a^4", w.deg(4) == 1)
chi = sum((-1) ** k for k in range(5))
C.check("consistency: w4[RP^4] = Euler characteristic mod 2 (= 1)", w.deg(4) == chi % 2, f"chi = {chi}")
C.done()
