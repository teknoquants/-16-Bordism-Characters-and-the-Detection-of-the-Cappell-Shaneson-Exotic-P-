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
"""Lemma 5.7 - the gap is measured from the ground multiplet."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.7 (Which energy gap is meaningful)",
            "For odd n, a J_n-symmetric H has a degenerate ground space; the gap above the ground multiplet "
            "is min{E - E0 : E in spec(H), E > E0}.")
for n in (1, 3, 5):
    d = 2 ** n
    U = kron(*([Y] * n))
    A = random_hermitian(d)
    Hm = (A + U @ np.conj(A) @ U.conj().T) / 2
    E = np.linalg.eigvalsh(Hm)
    naive = E[1] - E[0]
    distinct = np.unique(np.round(E, 9))
    C.check(f"n={n}: naive E1 - E0 = 0 (degenerate ground space)", abs(naive) < 1e-8)
    if len(distinct) == 1:
        C.check(f"n={n}: H is scalar (one Kramers pair), so no excited eigenvalue exists - the lemma's "
                "hypothesis 'if an excited eigenvalue exists' is genuinely needed", np.allclose(Hm, E[0] * np.eye(d)))
    else:
        gap = distinct[1] - distinct[0]
        C.check(f"n={n}: multiplet gap = min{{E - E0 : E > E0}} > 0", gap > 1e-6, f"gap = {gap:.4f}")
C.done()
