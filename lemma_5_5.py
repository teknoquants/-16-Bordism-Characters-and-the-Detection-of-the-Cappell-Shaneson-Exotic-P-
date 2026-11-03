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
"""Lemma 5.5 - [H] has order 2 in Br(R); Jn^2 = (-1)^n."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.5 (Brauer order and tensor parity)",
            "[H] has order two in Br(R), and Jn = Y^{⊗n} K satisfies Jn^2 = (-1)^n I.")
ok = True
for _ in range(500):
    p, q = RNG.normal(size=4), RNG.normal(size=4)
    ok &= np.isclose(np.linalg.norm(qmul(p, q)), np.linalg.norm(p) * np.linalg.norm(q))
C.check("H has no zero divisors (|pq| = |p||q|), so H is not split (H ≇ M2(R))", ok)
E12 = np.array([[0, 1], [0, 0]])
C.check("control: M2(R) has a non-zero nilpotent (zero divisor)", np.allclose(E12 @ E12, 0))
C.note("[H]^2 = [M4(R)] = [R] (Lemma 5.4), so the order is exactly 2", "from Lemma 5.4")
for n in range(1, 7):
    Yn = kron(*([Y] * n))
    Jn = Antilinear(Yn)
    C.check(f"n={n}: conj(Y^⊗n) = (-1)^n Y^⊗n and Jn^2 = (-1)^n I",
            np.allclose(np.conj(Yn), (-1) ** n * Yn) and np.allclose(Jn.square(), (-1) ** n * np.eye(2 ** n)))
# commutant of Jn: real dimension 4^n (M_{2^n}(R) for even n, M_{2^{n-1}}(H) for odd n)
for n in (1, 2, 3):
    d = 2 ** n
    Yn = kron(*([Y] * n))
    def vec(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
    def unvec(x): return (x[:d * d] + 1j * x[d * d:]).reshape(d, d)
    Lmat = np.column_stack([vec(unvec(e) - Yn @ np.conj(unvec(e)) @ np.linalg.inv(Yn)) for e in np.eye(2 * d * d)])
    C.check(f"n={n}: the real algebra commuting with Jn has dimension 4^n = {4 ** n}",
            2 * d * d - np.linalg.matrix_rank(Lmat) == 4 ** n)
C.done()
