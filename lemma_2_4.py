#!/usr/bin/env python3
"""Lemma 2.4 - Pin+ obstruction w2, Pin- obstruction w2 + w1^2."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 2.4 (Why the two Pin obstructions differ)",
            "E admits a Pin+ lift iff w2(E) = 0 and a Pin- lift iff w2(E) + w1(E)^2 = 0.",
            "Verifies the Clifford sign cocycle, its cocycle property, its identification with "
            "e2 + eps*e1^2 on (Z/2)^r, the reflection property of the e_i, and the symmetric-"
            "polynomial identities. The Čech/obstruction-theory step is standard and not simulated.")

for eps in (0, 1):
    sign = +1 if eps == 0 else -1
    for r in range(1, 6):
        e = clifford_generators(r, sign)
        dim = e[0].shape[0]
        okrel = all(np.allclose(e[i] @ e[j] + e[j] @ e[i], 2 * sign * np.eye(dim) * (i == j))
                    for i in range(r) for j in range(r))
        def s(g):
            M = np.eye(dim, dtype=complex)
            for i, gi in enumerate(bits(g, r)):
                if gi:
                    M = M @ e[i]
            return M
        def c_formula(g, h):
            gb, hb = bits(g, r), bits(h, r)
            return (sum(gb[i] * hb[j] for i in range(r) for j in range(i)) +
                    eps * sum(gb[i] * hb[i] for i in range(r))) % 2
        ok = True
        for g in range(1 << r):
            for h in range(1 << r):
                lhs, rhs = s(g) @ s(h), s(g ^ h)
                lam = (-1) ** c_formula(g, h)
                ok &= np.allclose(lhs, lam * rhs)
        name = "Pin+" if eps == 0 else "Pin-"
        C.check(f"{name}, r={r}: Clifford relations e_i e_j + e_j e_i = 2(±1)δ_ij", okrel)
        C.check(f"{name}, r={r}: s(g)s(h) = (-1)^c(g,h) s(g+h) with the stated c", ok)
        if r <= 4:
            C.check(f"{name}, r={r}: c is a 2-cocycle on (Z/2)^r", is_2cocycle(c_formula, r))
        # c is literally the cup-product cocycle of sum_{i>j} x_i x_j + eps sum x_i^2
        cup = lambda g, h: (sum(bits(g, r)[i] * bits(h, r)[j] for i in range(r) for j in range(i))
                            + eps * sum(bits(g, r)[i] * bits(h, r)[i] for i in range(r))) % 2
        C.check(f"{name}, r={r}: c equals the cup cocycle of e2(x) + eps*e1(x)^2",
                all(cup(g, h) == c_formula(g, h) for g in range(1 << r) for h in range(1 << r)))
        if 2 <= r <= 4:
            C.check(f"{name}, r={r}: the class is non-zero (c is not a coboundary)",
                    not is_2coboundary(c_formula, r))

# Reflection property: the twisted adjoint action of e_i reflects in e_i^perp (both signs).
for sign in (+1, -1):
    r = 4
    e = clifford_generators(r, sign)
    ok = True
    for i in range(r):
        v = RNG.normal(size=r)
        V = sum(v[j] * e[j] for j in range(r))
        W = -e[i] @ V @ np.linalg.inv(e[i])          # twisted adjoint alpha(e_i) V e_i^{-1}
        vref = v.copy(); vref[i] *= -1
        ok &= np.allclose(W, sum(vref[j] * e[j] for j in range(r)))
    C.check(f"e_i^2 = {sign:+d}: twisted adjoint of e_i is the reflection in e_i^⊥", ok)

# Symmetric-polynomial identities over F2: Π(1+x_i) has degree-2 part e2; (Σx_i)^2 = Σx_i^2.
def poly_mul(p, q):
    out = {}
    for m1 in p:
        for m2 in q:
            m = tuple(x + y for x, y in zip(m1, m2))
            out[m] = out.get(m, 0) ^ 1
    return {m for m, v in out.items() if v}
for r in range(2, 6):
    unit = tuple([0] * r)
    prod = {unit}
    for i in range(r):
        xi = tuple(1 if j == i else 0 for j in range(r))
        prod = poly_mul(prod, {unit, xi})
    e2 = {m for m in prod if sum(m) == 2}
    target = {tuple(1 if j in (i, k) else 0 for j in range(r)) for i in range(r) for k in range(i)}
    C.check(f"r={r}: degree-2 part of Π(1+x_i) is e2 = Σ_(i>j) x_i x_j", e2 == target)
    e1 = {tuple(1 if j == i else 0 for j in range(r)) for i in range(r)}
    sq = poly_mul(e1, e1)
    C.check(f"r={r}: (Σ x_i)^2 = Σ x_i^2 over F2",
            sq == {tuple(2 if j == i else 0 for j in range(r)) for i in range(r)})
C.done()
