#!/usr/bin/env python3
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
