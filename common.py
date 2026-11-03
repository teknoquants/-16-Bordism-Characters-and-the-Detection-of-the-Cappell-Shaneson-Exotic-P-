"""Shared helpers for the Problem 5 white-paper validation scripts.

Only the Python standard library and NumPy are used.
Every script imports `Checker`, records named checks, and exits with
status 1 if any check fails, so `run_all.py` can aggregate results.
"""
from __future__ import annotations

import itertools
import sys
from math import comb

import numpy as np

TOL = 1e-10
RNG = np.random.default_rng(20261103)


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------
class Checker:
    def __init__(self, label: str, statement: str, scope: str = ""):
        self.label = label
        self.passed = 0
        self.failed = 0
        print("=" * 78)
        print(label)
        print("-" * 78)
        print(statement.strip())
        if scope:
            print("\nScope: " + scope.strip())
        print("-" * 78)

    def check(self, name: str, cond: bool, detail: str = "") -> bool:
        cond = bool(cond)
        tag = "PASS" if cond else "FAIL"
        print(f"[{tag}] {name}" + (f"  ({detail})" if detail else ""))
        if cond:
            self.passed += 1
        else:
            self.failed += 1
        return cond

    def note(self, name: str, detail: str = "") -> None:
        """Record a step that is an external input or a structural remark (not counted as a check)."""
        print(f"[NOTE] {name}" + (f"  ({detail})" if detail else ""))

    def done(self) -> None:
        total = self.passed + self.failed
        print("-" * 78)
        print(f"{self.label}: {self.passed}/{total} checks passed")
        sys.exit(0 if self.failed == 0 else 1)


# --------------------------------------------------------------------------
# F2 arithmetic and truncated polynomial rings F2[a]/(a^{N+1})
# --------------------------------------------------------------------------
def binom2(n: int, k: int) -> int:
    """Binomial coefficient mod 2 (0 outside 0<=k<=n)."""
    if k < 0 or k > n or n < 0:
        return 0
    return comb(n, k) % 2


class TruncPoly:
    """Element of F2[a]/(a^{top+1}); coefficient list indexed by degree."""

    def __init__(self, coeffs, top: int):
        c = [int(x) % 2 for x in coeffs][: top + 1]
        c += [0] * (top + 1 - len(c))
        self.c, self.top = c, top

    @classmethod
    def gen(cls, top):
        return cls([0, 1], top)

    @classmethod
    def one(cls, top):
        return cls([1], top)

    def __add__(self, o):
        return TruncPoly([x ^ y for x, y in zip(self.c, o.c)], self.top)

    def __mul__(self, o):
        out = [0] * (self.top + 1)
        for i, x in enumerate(self.c):
            if x:
                for j, y in enumerate(o.c):
                    if y and i + j <= self.top:
                        out[i + j] ^= 1
        return TruncPoly(out, self.top)

    def __pow__(self, n):
        r = TruncPoly.one(self.top)
        for _ in range(n):
            r = r * self
        return r

    def __eq__(self, o):
        return self.c == o.c

    def deg(self, k):
        return self.c[k] if k <= self.top else 0

    def __repr__(self):
        terms = [("1" if k == 0 else ("a" if k == 1 else f"a^{k}"))
                 for k, x in enumerate(self.c) if x]
        return " + ".join(terms) if terms else "0"


def sq_power_of_a(m: int, k: int, top: int) -> TruncPoly:
    """Steenrod square Sq^k(a^m) = C(m,k) a^{m+k} in H*(RP^top; F2)."""
    c = [0] * (top + 1)
    if m + k <= top:
        c[m + k] = binom2(m, k)
    return TruncPoly(c, top)


# --------------------------------------------------------------------------
# Linear algebra over F2
# --------------------------------------------------------------------------
def f2_rank(M) -> int:
    A = (np.array(M, dtype=np.int64) % 2).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
        if r == rows:
            break
    return r


def f2_solvable(A, b) -> bool:
    """Is A x = b solvable over F2?"""
    A = np.array(A, dtype=np.int64) % 2
    b = np.array(b, dtype=np.int64).reshape(-1, 1) % 2
    return f2_rank(A) == f2_rank(np.hstack([A, b]))


# --------------------------------------------------------------------------
# Group cohomology of elementary abelian 2-groups G=(Z/2)^r (elements = bitmasks)
# --------------------------------------------------------------------------
def is_2cocycle(c, r: int) -> bool:
    n = 1 << r
    for g in range(n):
        for h in range(n):
            for k in range(n):
                # (delta c)(g,h,k) = c(h,k) + c(g+h,k) + c(g,h+k) + c(g,h)
                if (c(h, k) ^ c(g ^ h, k) ^ c(g, h ^ k) ^ c(g, h)) & 1:
                    return False
    return True


def is_2coboundary(c, r: int) -> bool:
    """Is c = delta f for some f: G -> F2, (delta f)(g,h) = f(g)+f(h)+f(g+h)?"""
    n = 1 << r
    A, b = [], []
    for g in range(n):
        for h in range(n):
            row = [0] * n
            row[g] ^= 1
            row[h] ^= 1
            row[g ^ h] ^= 1
            A.append(row)
            b.append(c(g, h) & 1)
    return f2_solvable(A, b)


def bits(g: int, r: int):
    return [(g >> i) & 1 for i in range(r)]


# --------------------------------------------------------------------------
# Matrices: Paulis, Clifford generators, antilinear operators
# --------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


def kron(*ms):
    out = np.array([[1]], dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


def pauli_string(s: str):
    return kron(*[PAULI[ch] for ch in s])


def clifford_generators(r: int, sign: int):
    """r pairwise anticommuting matrices with e_i^2 = sign*I (Jordan-Wigner)."""
    m = (r + 1) // 2
    gens = []
    for k in range(m):
        for P in (X, Y):
            gens.append(kron(*([Z] * k + [P] + [I2] * (m - k - 1))))
    gens = gens[:r]
    if sign == -1:
        gens = [1j * g for g in gens]
    return gens


class Antilinear:
    """The antilinear operator v -> U conj(v)."""

    def __init__(self, U):
        self.U = np.array(U, dtype=complex)

    def __call__(self, v):
        return self.U @ np.conj(v)

    def __matmul__(self, other):
        if isinstance(other, Antilinear):          # (U K)(V K) = U conj(V)  (linear)
            return self.U @ np.conj(other.U)
        return Antilinear(self.U @ np.conj(other))  # (U K) V = U conj(V) K

    def square(self):
        return self.U @ np.conj(self.U)


def random_hermitian(d: int):
    A = RNG.normal(size=(d, d)) + 1j * RNG.normal(size=(d, d))
    return (A + A.conj().T) / 2


def random_complex(shape):
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)


# --------------------------------------------------------------------------
# Quaternions as numpy 4-vectors (w, x, y, z)
# --------------------------------------------------------------------------
def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([
        a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
        a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
        a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
        a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    ])


def qconj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])


QBASIS = [np.eye(4)[i] for i in range(4)]  # 1, i, j, k


def roots16(k: int):
    return np.exp(2j * np.pi * k / 16)


__all__ = [n for n in dir() if not n.startswith("_")] + ["itertools", "np", "comb"]
