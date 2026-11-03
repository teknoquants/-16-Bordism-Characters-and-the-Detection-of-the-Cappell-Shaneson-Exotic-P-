#!/usr/bin/env python3
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
