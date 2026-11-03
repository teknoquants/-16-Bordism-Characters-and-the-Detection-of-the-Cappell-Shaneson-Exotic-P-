#!/usr/bin/env python3
"""Lemma A.2 - T^2 = P gives real even sector and quaternionic odd sector."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma A.2 (The symmetry motivation is compatible at the level of signs)",
            "If T is antilinear with T^2 = P, then T preserves the grading, T^2 = I on the even sector and "
            "T^2 = -I on the odd sector, where i and T generate a quaternion action.")
de, do = 2, 4
d = de + do
P = np.diag([1] * de + [-1] * do).astype(complex)
U0 = np.zeros((d, d), complex)
U0[:de, :de] = np.eye(de)
U0[de:, de:] = kron(I2, Y)
# make it generic: conjugate by a random grading-preserving unitary W:  T = W T0 W^{-1} = W U0 conj(W)^{-1} K
def rand_unitary(k):
    q, _ = np.linalg.qr(random_complex((k, k))); return q
W = np.zeros((d, d), complex); W[:de, :de] = rand_unitary(de); W[de:, de:] = rand_unitary(do)
T = Antilinear(W @ U0 @ np.linalg.inv(np.conj(W)))
T2 = T.square()
C.check("T^2 = P", np.allclose(T2, P))
TP = T @ P; PT = Antilinear(P @ T.U)
C.check("TP = PT", np.allclose(TP.U, PT.U))
ve = np.concatenate([random_complex(de), np.zeros(do)])
vo = np.concatenate([np.zeros(de), random_complex(do)])
C.check("T maps even vectors to even vectors", np.allclose(T(ve)[de:], 0))
C.check("T maps odd vectors to odd vectors", np.allclose(T(vo)[:de], 0))
C.check("T^2 = +I on the even sector", np.allclose(T2[:de, :de], np.eye(de)))
C.check("T^2 = -I on the odd sector", np.allclose(T2[de:, de:], -np.eye(do)))
To = Antilinear(T.U[de:, de:])
iT = Antilinear(1j * To.U)
v = random_complex(do)
C.check("odd sector: T i = -i T", np.allclose(To(1j * v), -1j * To(v)))
C.check("odd sector: (iT)^2 = -I", np.allclose(iT.square(), -np.eye(do)))
C.check("odd sector: with I = i, J = T, K = iT:  J·K = I  (T(iT) = -i T^2 = i)", np.allclose(To @ iT, 1j * np.eye(do)))
C.check("odd sector: K·I = J  ((iT) i = T)", np.allclose((iT @ (1j * np.eye(do))).U, To.U))
C.check("odd sector: I^2 = J^2 = K^2 = -1", np.allclose(To.square(), -np.eye(do)) and np.allclose(iT.square(), -np.eye(do)))
C.check("odd-sector dimension is even (quaternionic)", do % 2 == 0)
C.done()
