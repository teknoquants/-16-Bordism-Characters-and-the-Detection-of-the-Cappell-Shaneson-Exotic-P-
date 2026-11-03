#!/usr/bin/env python3
"""Lemma 2.1 - H*(RP^4;F2) = F2[a]/(a^5), a = w1(gamma)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 2.1 (Cohomology and the tautological class)",
            "H*(RP^4;F2) = F2[a]/(a^5), deg a = 1, a = w1(gamma); a and a^2 are nonzero.")

# (1) Cellular cochains of RP^n: delta^k : C^k -> C^{k+1} is multiplication by 1+(-1)^{k+1}.
n = 4
deltas = [1 + (-1) ** (k + 1) for k in range(n)]           # k = 0..3
C.check("integral cellular coboundaries are 0 or 2", all(d in (0, 2) for d in deltas), str(deltas))
C.check("all coboundaries vanish mod 2", all(d % 2 == 0 for d in deltas))
dims = []
for k in range(n + 1):
    ker = 1 if (k == n or deltas[k] % 2 == 0) else 0
    im = 1 if (k > 0 and deltas[k - 1] % 2 == 1) else 0
    dims.append(ker - im)
C.check("dim H^k(RP^4;F2) = 1 for k = 0..4", dims == [1] * 5, str(dims))

# (2) Normalized bar complex of G = Z/2 = {1, g} with F2 coefficients.
def bar_boundary(chain):
    """Boundary of [g1|...|gk] (entries in {0,1}, 1 = g) in the normalized bar complex mod 2.
    Returns dict of chains with coefficients mod 2; degenerate chains (with an identity) vanish."""
    k = len(chain)
    out = {}
    def add(c):
        if 0 in c:            # normalized: chains containing the identity are zero
            return
        out[c] = out.get(c, 0) ^ 1
    add(tuple(chain[1:]))                                  # g1 . [g2|...|gk], trivial action
    for i in range(k - 1):
        add(tuple(chain[:i]) + ((chain[i] + chain[i + 1]) % 2,) + tuple(chain[i + 2:]))
    add(tuple(chain[:-1]))                                 # (-1)^k [g1|...|g_{k-1}]
    return {c: v for c, v in out.items() if v}

for k in range(1, 8):
    C.check(f"[g|...|g] (k={k}) is a cycle mod 2", bar_boundary([1] * k) == {})

def cup_power_eval(k):
    """Alexander-Whitney: (u ∪ ... ∪ u)[g1|...|gk] = u(g1)...u(gk), u(g) = 1."""
    u = {0: 0, 1: 1}
    val = 1
    for x in [1] * k:
        val *= u[x]
    return val % 2
for k in range(1, 8):
    C.check(f"u^{k} evaluates to 1 on the cycle [g|...|g], so u^{k} is not a coboundary",
            cup_power_eval(k) == 1)
C.check("the normalized bar complex has exactly one nondegenerate chain per degree",
        all(len([c for c in itertools.product([0, 1], repeat=k) if 0 not in c]) == 1 for k in range(1, 8)))

# (3) Restriction to the 4-skeleton: truncated ring.
a = TruncPoly.gen(4)
C.check("a^4 != 0 and a^5 = 0 in F2[a]/(a^5)", (a ** 4).c[4] == 1 and (a ** 5) == TruncPoly([0], 4))
C.check("a != 0 and a^2 != 0", (a ** 1).c[1] == 1 and (a ** 2).c[2] == 1)

# (4) w1(gamma) != 0: the tautological line over RP^1 is the Moebius band.
theta = np.linspace(0, 2 * np.pi, 2001)
frame = np.stack([np.cos(theta / 2), np.sin(theta / 2)], axis=1)   # spans the line at angle theta/2
continuous = np.all(np.linalg.norm(np.diff(frame, axis=0), axis=1) < 1e-2)
C.check("a continuous frame along the generator loop returns with sign -1 (non-trivial monodromy)",
        continuous and np.allclose(frame[-1], -frame[0]))
C.done()
