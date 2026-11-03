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
"""Lemma 5.3 - Phi(a⊗b)(x) = a x conj(b) is a unital algebra homomorphism."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.3 (A real algebra homomorphism)",
            "Φ(a⊗b)(x) = a x b̄ defines a unital real algebra homomorphism H ⊗_R H -> End_R(H).")
def Phi(a, b):
    return np.column_stack([qmul(qmul(a, e), qconj(b)) for e in QBASIS])
C.check("Hamilton relations: i^2 = j^2 = k^2 = ijk = -1",
        all(np.allclose(qmul(q, q), -QBASIS[0]) for q in QBASIS[1:]) and
        np.allclose(qmul(qmul(QBASIS[1], QBASIS[2]), QBASIS[3]), -QBASIS[0]))
C.check("Φ(1⊗1) = id", np.allclose(Phi(QBASIS[0], QBASIS[0]), np.eye(4)))
ok = all(np.allclose(Phi(qmul(a, c), qmul(b, d)), Phi(a, b) @ Phi(c, d))
         for a in QBASIS for b in QBASIS for c in QBASIS for d in QBASIS)
C.check("Φ(ac⊗bd) = Φ(a⊗b)Φ(c⊗d) on all 256 basis quadruples", ok)
ok = True
for _ in range(200):
    a, b, c, d = (RNG.normal(size=4) for _ in range(4))
    ok &= np.allclose(Phi(qmul(a, c), qmul(b, d)), Phi(a, b) @ Phi(c, d))
C.check("multiplicativity on 200 random quadruples", ok)
C.check("conj(bd) = conj(d) conj(b)", all(np.allclose(qconj(qmul(b, d)), qmul(qconj(d), qconj(b)))
                                          for b in QBASIS for d in QBASIS))
a, b, x = RNG.normal(size=4), RNG.normal(size=4), RNG.normal(size=4)
C.check("real bilinearity in a and b", np.allclose(Phi(2 * a + b, b), 2 * Phi(a, b) + Phi(b, b)))
C.done()
