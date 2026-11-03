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
"""Lemma 5.2 - matrices commuting with J1 are [[alpha,-conj(beta)],[beta,conj(alpha)]]."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.2 (Fixed matrices)",
            "The complex 2x2 matrices commuting with J1 are exactly [[α, -β̄],[β, ᾱ]].")
# Real-linear map L(M) = M - Y conj(M) Y on M2(C) ≅ R^8; its kernel is the commutant.
def vec(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
def unvec(x): return (x[:4] + 1j * x[4:]).reshape(2, 2)
L = np.column_stack([vec(unvec(e) - Y @ np.conj(unvec(e)) @ Y) for e in np.eye(8)])
null_dim = 8 - np.linalg.matrix_rank(L)
C.check("the commutant has real dimension 4 (= dim H)", null_dim == 4)
a, b, c, d = sp = [complex(*RNG.normal(size=2)) for _ in range(4)]
M = np.array([[a, b], [c, d]])
C.check("Y conj(M) Y = [[d̄, -c̄],[-b̄, ā]]",
        np.allclose(Y @ np.conj(M) @ Y, np.array([[np.conj(d), -np.conj(c)], [-np.conj(b), np.conj(a)]])))
ok = True
for _ in range(300):
    al, be = random_complex(2)
    Q = np.array([[al, -np.conj(be)], [be, np.conj(al)]])
    ok &= np.allclose(Y @ np.conj(Q) @ Y, Q)
C.check("every matrix of the stated form commutes with J1", ok)
basis = [np.array([[1, 0], [0, 1]]), np.array([[1j, 0], [0, -1j]]),
         np.array([[0, -1], [1, 0]]), np.array([[0, 1j], [1j, 0]])]
C.check("the stated form spans the full 4-dim commutant",
        np.linalg.matrix_rank(np.column_stack([vec(B.astype(complex)) for B in basis])) == 4
        and all(np.allclose(L @ vec(B.astype(complex)), 0) for B in basis))
M2 = random_complex((2, 2))
C.check("control: a random matrix does not commute with J1", not np.allclose(Y @ np.conj(M2) @ Y, M2))
C.done()
