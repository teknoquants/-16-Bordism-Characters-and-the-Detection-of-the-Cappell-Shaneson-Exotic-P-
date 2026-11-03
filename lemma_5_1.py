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
"""Lemma 5.1 - J1^2 = -I and J1 M J1^-1 = Y conj(M) Y."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.1 (Square and conjugation formula)",
            "J1^2 = -I, and J1 M J1^{-1} = Y conj(M) Y for every complex 2x2 matrix M.")
J1 = Antilinear(Y)
C.check("conj(Y) = -Y and Y^2 = I", np.allclose(np.conj(Y), -Y) and np.allclose(Y @ Y, I2))
C.check("J1^2 = Y conj(Y) = -I", np.allclose(J1.square(), -I2))
J1inv = Antilinear(-Y)                       # (YK)^{-1} = KY = conj(Y) K = -Y K
C.check("(YK)^{-1} = KY = -YK", np.allclose(J1 @ J1inv, I2) and np.allclose(J1inv @ J1, I2))
ok = True
for _ in range(500):
    M = random_complex((2, 2))
    conjugated = J1 @ Antilinear(M @ (-Y))   # J1 ∘ M ∘ J1^{-1}  (M ∘ (-Y K) = (-M Y) K)
    ok &= np.allclose(conjugated, Y @ np.conj(M) @ Y)
C.check("J1 M J1^{-1} = Y conj(M) Y for 500 random M", ok)
v = random_complex(2)
C.check("antilinearity: J1(i v) = -i J1(v)", np.allclose(J1(1j * v), -1j * J1(v)))
C.done()
