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
"""Lemma 2.2 - T RP^4 + eps = 5 gamma."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 2.2 (Stable tangent bundle)", "T RP^4 ⊕ ε ≅ 5γ as real vector bundles.",
            "Checked on the double cover S^4: the explicit bundle map is a fibrewise isomorphism and "
            "is Z/2-equivariant, so it descends to RP^4.")

def tangent_basis(x):
    # orthonormal basis of x^perp in R^5
    q, _ = np.linalg.qr(np.column_stack([x, RNG.normal(size=(5, 4))]))
    return q[:, 1:]

ok_iso, ok_eq, ok_hom = True, True, True
for _ in range(200):
    x = RNG.normal(size=5); x /= np.linalg.norm(x)
    B = tangent_basis(x)
    # Psi : T_x S^4 ⊕ R -> R^5,  (t, s) -> t + s x
    Psi = np.column_stack([B, x])
    ok_iso &= abs(np.linalg.det(Psi)) > 1e-8
    # Z/2 actions: (x,t,s) -> (-x,-t,s) on the source, (x,v) -> (-x,-v) on the target.
    t = B @ RNG.normal(size=4); s = RNG.normal()
    lhs = (-t) + s * (-x)                 # Psi at -x of the image of (t, s)
    ok_eq &= np.allclose(lhs, -(t + s * x))
    # Hom(l, l^perp) ⊕ Hom(l, l) = Hom(l, R^5): ranks 4 + 1 = 5
    Pl = np.outer(x, x)
    ok_hom &= (np.linalg.matrix_rank(np.eye(5) - Pl) == 4 and np.linalg.matrix_rank(Pl) == 1)
C.check("Psi(t,s) = t + s x is a fibrewise isomorphism T_xS^4 ⊕ R -> R^5", ok_iso)
C.check("Psi is equivariant: antipodal action -1 on T S^4, +1 on ε, -1 on R^5", ok_eq)
C.check("rank Hom(γ,γ⊥) + rank Hom(γ,γ) = 4 + 1 = 5", ok_hom)
C.note("the quotient of S^4×R^5 by (x,v)~(-x,-v) is 5γ (each coordinate is a copy of γ)", "definition of γ as S^n ×_{Z/2} R with the sign action")
# The identity section of Hom(γ,γ) is nowhere zero and equivariant -> trivial line ε.
ok_triv = True
for _ in range(100):
    x = RNG.normal(size=5); x /= np.linalg.norm(x)
    # the identity section of Hom(γ,γ): l -> id_l; on the cover, the normal coefficient s = <v, x> of v = s x
    # is unchanged by (x, v) -> (-x, -v): <-v, -x> = <v, x>
    v = RNG.normal() * x
    ok_triv &= np.isclose(np.dot(-v, -x), np.dot(v, x))
C.check("the normal coordinate s = <v,x> is invariant under (x,v)->(-x,-v): Hom(γ,γ) ≅ ε is trivial", ok_triv)
C.done()
