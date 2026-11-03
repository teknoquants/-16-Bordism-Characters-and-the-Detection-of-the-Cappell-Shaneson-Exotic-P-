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
"""Lemma 5.6 - Kramers degeneracy."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 5.6 (Kramers degeneracy)",
            "If H = H^* commutes with an antiunitary J with J^2 = -I, every eigenspace has even dimension.")
for n in (1, 3, 5):
    d = 2 ** n
    U = kron(*([Y] * n))
    J = Antilinear(U)
    C.check(f"n={n}: J^2 = -I", np.allclose(J.square(), -np.eye(d)))
    A = random_hermitian(d)
    Hm = (A + U @ np.conj(A) @ U.conj().T) / 2          # J-symmetrised Hamiltonian
    C.check(f"n={n}: H is Hermitian and commutes with J",
            np.allclose(Hm, Hm.conj().T) and np.allclose(U @ np.conj(Hm), Hm @ U))
    E, V = np.linalg.eigh(Hm)
    # group eigenvalues
    groups = []
    for e in E:
        if groups and abs(e - groups[-1][0]) < 1e-8:
            groups[-1][1] += 1
        else:
            groups.append([e, 1])
    C.check(f"n={n}: every eigenvalue has even multiplicity", all(m % 2 == 0 for _, m in groups),
            f"multiplicities {[m for _, m in groups][:6]}...")
    v = V[:, 0]
    Jv = J(v)
    C.check(f"n={n}: <Jv, v> = 0 and HJv = E Jv for a ground vector",
            abs(np.vdot(Jv, v)) < 1e-9 and np.allclose(Hm @ Jv, E[0] * Jv))
    x, y = random_complex(d), random_complex(d)
    C.check(f"n={n}: antiunitarity <Jx,Jy> = <y,x>", np.isclose(np.vdot(J(x), J(y)), np.vdot(y, x)))
# control: with J^2 = +I (n = 2) degeneracy is not forced
U2 = kron(Y, Y); A = random_hermitian(4); H2 = (A + U2 @ np.conj(A) @ U2.conj().T) / 2
E2 = np.linalg.eigvalsh(H2)
C.check("control n=2 (J^2 = +I): spectrum is generically non-degenerate", np.min(np.diff(E2)) > 1e-6)
C.done()
