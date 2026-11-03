# Problem 5 White Paper — Validation Scripts

This folder holds the validation scripts for *Problem 5: Pin⁺ Exotic Detection*, corrected edition (second review), 3 November 2026. There is one script for every lemma and the theorem: 31 lemmas plus Theorem 3.6, so 32 scripts.

## Running

```bash
pip install numpy          # the only dependency (Python 3.9+)
python run_all.py          # summary table; exit status 0 only if all scripts pass
python run_all.py -v       # full output of every check
python scripts/lemma_2_4.py   # run a single script
```

Each script prints the lemma's statement and its scope. It then runs a list of named `[PASS]`/`[FAIL]` checks and exits with status 1 if any check fails. `[NOTE]` lines record external inputs or structural remarks and are not counted as checks. `RESULTS.txt` and `RESULTS_verbose.txt` hold the output of the reference run. Random tests use a fixed seed (20261103), so runs are reproducible.

## What each script validates

| Script | Statement | Method |
|---|---|---|
| `lemma_1_1` | T² = P is not an involution | explicit antilinear T on a super space |
| `lemma_2_1` | H\*(ℝP⁴;𝔽₂) = 𝔽₂[a]/(a⁵) | cellular coboundaries; normalized bar complex of ℤ/2 (cycles, cup powers); Möbius monodromy |
| `lemma_2_2` | TℝP⁴ ⊕ ε ≅ 5γ | explicit bundle map on S⁴, fibrewise isomorphism and ℤ/2-equivariance (random points) |
| `lemma_2_3` | w = 1 + a + a⁴ | truncated polynomial arithmetic over 𝔽₂; w₄ = χ mod 2 |
| `lemma_2_4` | Pin⁺ ⇔ w₂ = 0; Pin⁻ ⇔ w₂ + w₁² = 0 | explicit Clifford matrices (r ≤ 5, both signs); sign cocycle vs formula; cocycle and non-coboundary tests over 𝔽₂; reflection property; symmetric-polynomial identities |
| `lemma_2_5` | non-orientable, Pin⁺, no Pin⁻ | two routes: characteristic classes, and lifting the holonomy −I₅ into Pin±(5) |
| `lemma_2_6` | exactly two Pin⁺ structures | H¹-torsor count; the two lifts ±e₁⋯e₅ |
| `lemma_2_7` | same for any homotopy ℝP⁴ | SW classes from Wu's formula (homotopy-invariant data), cross-checked against (1+a)ⁿ⁺¹ for ℝPⁿ, n ≤ 10 |
| `lemma_3_1` | χ(k) = ζᵏ injective character | exhaustive over ℤ/16 |
| `lemma_3_2` | I = χ([M,p]) invariant, multiplicative | bordism modelled as ℤ/16 |
| `lemma_3_3` | separation for all labellings | all units of ℤ/16, all injective characters; control: even u fail |
| `lemma_3_4` | A(R) = cos(π/8) = −A(Q) | all odd u; non-vanishing |
| `lemma_3_5` | functorial extension | algebraic layer only: ω = χ∘s is a character for random abelian models of s |
| `theorem_3_6` | Pin⁺ detection | assembles 2.6, 2.7, 3.1, 3.3; inputs flagged |
| `lemma_4_1` | oriented obstruction not applicable | w₁ ≠ 0 for R and Q |
| `lemma_4_2`–`4_4` | scope lemmas | illustrative counterexamples (16 distinct theories share the same data; conjugation gives χ⁻¹) |
| `lemma_4_5` | no detection of homotopy 4-spheres | Wu classes, signature, Rokhlin consistency (K3: σ = −16, Â = 2) |
| `lemma_5_1`–`5_2` | J₁² = −I; commutant of J₁ | random matrices; real-linear null space of dimension 4 |
| `lemma_5_3`–`5_4` | ℍ⊗ℍ ≅ M₄(ℝ) | all 256 basis quadruples; trace identity; Gram matrix = 4I; rank 16 |
| `lemma_5_5` | Brauer order 2; Jₙ² = (−1)ⁿ | norm multiplicativity; n = 1..6; commutant of Jₙ has real dimension 4ⁿ |
| `lemma_5_6`–`5_7` | Kramers degeneracy; multiplet gap | random Jₙ-symmetric Hamiltonians, n = 1, 3, 5; control at n = 2 |
| `lemma_A_1`–`A_5` | appendix lemmas | Pin± admissibility for R and Q; quaternion relations on the odd sector; a finite invariant list (χ, all SW numbers, Dijkgraaf–Witten counts for ℤ/2, ℤ/3, ℤ/4, S₃) agrees on R and Q while the Pin⁺ bordism character does not |
| `lemma_B_1` | arithmetic checks | direct |

## What is *not* validated

The scripts check every computation and every finite or algebraic step in the paper. They do not, and cannot, prove the paper's external inputs, which are taken from the literature:

- **Input B.** Ω₄^{Pin⁺} ≅ ℤ/16; ℝP⁴ ↦ ±1 and the Cappell–Shaneson ℝP⁴ ↦ ±9 (Stolz; Debray, Ex. 3.4; Kirby–Taylor); and that the two manifolds are homeomorphic (Hambleton–Kreck–Teichner). These values are hard-coded as `R_CLASSES` and `Q_CLASSES` and labelled in every script that uses them.
- **Input E.** Reutter–Schommer-Pries, Theorem 3.14: invertible super TQFTs ↔ characters of the semistable bordism group.
- **Rokhlin's theorem.** Ω₄^{Spin} ≅ ℤ, detected by σ/16 (used in Lemma 4.5).
- **Obstruction theory.** The Čech identification in Lemma 2.4 is standard; only its algebraic core (the extension class on (ℤ/2)ʳ) is computed.

Lemmas 4.2, 4.3, 4.4 and A.4 are scope statements. Their scripts give illustrative witnesses (explicit examples), not proofs. Passing scripts are evidence for the correctness of the computations, not a formal proof; a formalization (for example in Lean) would be the next step.

## Mutation check

To confirm the checks are not vacuous, three errors were planted and each was caught (non-zero exit):

- Cappell–Shaneson classes changed to ±1 → `lemma_3_3` fails.
- Pin± sign swapped in the Clifford cocycle → `lemma_2_4` fails.
- Wu class formula perturbed → `lemma_2_7` fails.

## Licence

MIT, consistent with the Challenge Series validation code.
