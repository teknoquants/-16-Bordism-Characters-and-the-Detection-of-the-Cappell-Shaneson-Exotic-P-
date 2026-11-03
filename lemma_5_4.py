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
"""Lemma 5.4 - Phi is an isomorphism H⊗H ≅ M4(R)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.4 (Injectivity and surjectivity)", "Φ is an isomorphism; H ⊗_R H ≅ M4(R).")
L = lambda u: np.column_stack([qmul(u, e) for e in QBASIS])
Rm = lambda v: np.column_stack([qmul(e, v) for e in QBASIS])
ok = all(np.isclose(np.trace(L(u) @ Rm(v)), 4 * u[0] * v[0]) for u in QBASIS for v in QBASIS)
C.check("tr(L_u R_v) = 4 Re(u) Re(v) for all basis pairs", ok)
ok = all(np.allclose(L(a).T, L(qconj(a))) and np.allclose(Rm(a).T, Rm(qconj(a))) for a in QBASIS)
C.check("adjoints: L_a^* = L_ā, R_a^* = R_ā", ok)
Phi = lambda a, b: L(a) @ Rm(qconj(b))
imgs = [Phi(a, b) for a in QBASIS for b in QBASIS]
G = np.array([[np.trace(A.T @ B) for B in imgs] for A in imgs])
C.check("Gram matrix tr(Φ(a⊗b)^*Φ(c⊗d)) = 4 δ_ac δ_bd", np.allclose(G, 4 * np.eye(16)))
C.check("the 16 images are linearly independent (rank 16)",
        np.linalg.matrix_rank(np.array([A.ravel() for A in imgs])) == 16)
C.check("dim_R(H⊗H) = 16 = dim_R M4(R)", 4 * 4 == 16)
C.done()
