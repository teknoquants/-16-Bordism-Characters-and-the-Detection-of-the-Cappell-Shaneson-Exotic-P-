#!/usr/bin/env python3
"""Lemma A.1 - the exercise selects the positive Pin sign."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma A.1 (The exercise selects the positive Pin sign)",
            "RP^4 (and, by Lemma 2.7, the Cappell-Shaneson RP^4) can be evaluated by a Pin+ theory and "
            "not by a Pin--only theory. Among Pin± theories, the test must use Pin+.")
a = TruncPoly.gen(4)
w = (TruncPoly.one(4) + a) ** 5          # same for any homotopy RP^4 (Lemma 2.7)
w1, w2 = w.deg(1), w.deg(2)
pin_plus = (w2 == 0)
pin_minus = ((w2 + w1 * w1) % 2 == 0)
for name in ("RP^4", "Cappell-Shaneson RP^4"):
    C.check(f"{name}: admissible for Pin+", pin_plus)
    C.check(f"{name}: not admissible for Pin-", not pin_minus)
C.check("among {Pin+, Pin-}, only Pin+ admits both manifolds",
        [s for s, ok in (("Pin+", pin_plus), ("Pin-", pin_minus)) if ok] == ["Pin+"])
C.done()
