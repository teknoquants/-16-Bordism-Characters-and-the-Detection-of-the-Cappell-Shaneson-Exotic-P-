#!/usr/bin/env python3
"""Lemma 4.5 - the theory of Theorem 3.6 cannot see homotopy 4-spheres."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import *  # noqa: F401,F403

C = Checker("Lemma 4.5 (The theory of Theorem 3.6 cannot see homotopy 4-spheres)",
            "A smooth homotopy 4-sphere Σ has a unique Pin+ structure and Z_ω(Σ) = 1 = Z_ω(S^4).",
            "Uses Rokhlin's theorem / Ω4^Spin ≅ Z via σ/16 (external). Checks the cohomological and "
            "arithmetic steps.")
# H*(Σ;F2): F2 in degrees 0 and 4 only.
H = [1, 0, 0, 0, 1]
# Wu classes: v_k with Sq^k: H^{4-k} -> H^4. H^3 = H^2 = 0 so v1 = v2 = 0, hence w = Sq(v) = 1.
v1 = 0 if H[3] == 0 else None
v2 = 0 if H[2] == 0 else None
C.check("Wu classes v1 = v2 = 0 (H^3 = H^2 = 0)", v1 == 0 and v2 == 0)
C.check("w1(Σ) = w2(Σ) = 0 => orientable, spin and Pin+", H[1] == 0 and H[2] == 0)
C.check("number of Pin+ structures = |H^1(Σ;F2)| = 1", 2 ** H[1] == 1)
sigma = 0  # H^2(Σ;Z) = 0 so the intersection form is empty
C.check("signature of Σ is 0 (empty intersection form)", sigma == 0)
# Spin bordism detected by signature/16 (Rokhlin). K3 consistency: σ = -16, p1 = 3σ, Â = -p1/24 = 2.
sigK3 = -16; p1 = 3 * sigK3; Ahat = -p1 / 24
C.check("Rokhlin consistency: K3 has σ = -16, Â = -p1/24 = 2 (integer, even)", Ahat == 2 and sigK3 % 16 == 0)
spin_class = sigma // 16
C.check("[Σ] = σ/16 · [K3] = 0 in Ω4^Spin, hence 0 in Ω4^Pin+", spin_class == 0)
C.check("Z_ω(Σ) = χ(0) = 1 = Z_ω(S^4)", np.isclose(roots16(0), 1))
C.done()
