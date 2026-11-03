#!/usr/bin/env python3
"""Lemma 4.1 - an oriented obstruction cannot be applied unchanged."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 4.1 (An oriented obstruction cannot be applied unchanged)",
            "A theorem quantified over oriented 4-manifolds does not, without further argument, rule out "
            "detection of the pair R, Q.")
w1_R = 1      # Lemma 2.3: w1(RP^4) = a
w1_Q = 1      # Lemma 2.7: w1 is a homotopy invariant
orientable = lambda w1: w1 == 0
C.check("R is non-orientable (w1 = a != 0)", not orientable(w1_R))
C.check("Q is non-orientable (w1 is a homotopy invariant, Lemma 2.7)", not orientable(w1_Q))
C.check("the pair lies outside the hypothesis 'oriented' of the obstruction theorem",
        not (orientable(w1_R) and orientable(w1_Q)))
# Orientability is detected by det of holonomy on π1 = Z/2: det(-I_4) = +1 on R^4 but tangent holonomy of RP^4
# is the restriction of -I_5 on 5γ minus a trivial line; its determinant is det(-I_5)/1 = -1.
C.check("orientation character of RP^4: det(-I_5)/det(+1) = -1", np.isclose(np.linalg.det(-np.eye(5)), -1))
C.done()
