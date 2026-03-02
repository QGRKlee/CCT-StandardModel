# Machine Verification Suite
## Kinematic Spinors on Division-Algebra Root Systems

**Paper:** "Kinematic Spinors on Division-Algebra Root Systems:
A First-Principles Derivation of the Standard Model Weyl Group"
— Klee Irwin (2026)

This repository contains machine-verifiable proofs for all five theorems
in the paper. Every computational claim prints `[PASS]` or `[FAIL]`;
the suite exits 0 only if all claims pass.

## Quick start

```bash
# Python 3.10+, only dependency is numpy
pip install -r requirements.txt

# Run all theorems (A, C, D, E finish in ~20 s; B takes ~10-15 min)
python run_all.py

# Quick check — skip Theorem B
python run_all.py --skip-partition

# Single theorem
python run_all.py --theorem A
```

## Repository layout

```
├── e8_utils.py             # Shared: root systems, Hopf map, permutation algebra
├── verify_hierarchy.py     # Theorem A — A₁⊂A₂⊂D₄⊂E₈ hierarchy + kinematic spinors
├── verify_partition.py     # Theorem B — 15,120 orbits, stabilizer, kernel ≅ 2T
├── verify_gut_breaking.py  # Theorem C — W(D₅)⊃W(A₄)⊃W(A₂)×W(A₁)
├── verify_triality.py      # Theorem D — D₄ triality, 8 partitions, S₃ action
├── verify_angles.py        # Theorem E — semi-conformal projection, angle preservation
├── run_all.py              # Master runner with summary table
├── requirements.txt        # numpy >= 1.24
└── README.md
```

## What each theorem verifies

### Theorem A — Root-system hierarchy (14 claims, <1 s)
- A₁, A₂, D₄, E₈ root systems: correct counts, norms, ranks
- Level-1 kinematic spinor: Spin(2) rotor of order 6, R³ = −I (spinorial signature)
- Compositional structure: A₂ = 3 × A₁, D₄ = 4 × A₂, E₈ = 10 × D₄ via Hopf

### Theorem B — Partition orbit and stabilizer (16 claims, ~10-15 min)
- 10 Hopf shells, 5 perpendicular pairs
- BFS orbit under W(E₈) reflections: exactly 15,120 partitions
- Orbit-stabilizer: |Stab| = 696,729,600 / 15,120 = 46,080
- Schreier generator closure reaches 46,080
- Kernel: 24 elements, order distribution {1:1, 2:1, 3:8, 4:6, 6:8} → 2T = SL(2,F₃)
- Image: 1,920 elements preserving perpendicular pairs → W(D₅)
- Extension 1 → 2T → Stab → W(D₅) → 1 is non-split

### Theorem C — GUT breaking chain (9 claims, <10 s)
- W(D₅) = (ℤ/2)⁴ ⋊ S₅ of order 1,920
- W(A₄) = S₅ subgroup of order 120, index 16
- W(A₂) × W(A₁) = S₃ × S₂ of order 12, index 10
- All 10 = C(5,3) Standard Model embeddings verified
- Even parity: all W(D₅) elements have det = +1

### Theorem D — D₄ triality (9+ claims, <1 s)
- 16 A₂ sub-root-systems of D₄ (exhaustive search)
- 8 partitions of D₄ into 4 disjoint A₂'s
- Triality generators τ (order 3) and σ (order 2) generate S₃
- S₃ acts on the 8 partitions: 2 orbits of sizes {2, 6}
- Burnside verification: Σ|Fix(g)|/6 = 2

### Theorem E — Semi-conformal angle preservation (8 claims, <5 s)
- D₄ angle set in 4D = {60°, 90°, 120°, 180°}
- Projection along Hopf fiber normal → two concentric shells
- Outer (cuboctahedral, 12 pts): angles ⊆ {60°, 90°, 120°, 180°}
- Inner (octahedral): angles ⊆ {90°, 180°}
- Result independent of projection direction within normal space

## Dependencies

Only **NumPy** (≥ 1.24). No other packages required.
Tested with Python 3.10–3.12 on Windows and Linux.

## Expected output

```
==================================================================
  MACHINE VERIFICATION SUITE
  Kinematic Spinors on Division-Algebra Root Systems
==================================================================

Running verify_hierarchy.py ...
  [PASS] ...
  ...

==================================================================
  VERIFICATION SUMMARY
==================================================================
  Theorem A (Hierarchy)            : 14/14 PASS   (0.3s)
  Theorem B (Partition)            : 16/16 PASS   (612.0s)
  Theorem C (GUT Breaking)         :  9/ 9 PASS   (3.2s)
  Theorem D (Triality)             :  9/ 9 PASS   (0.4s)
  Theorem E (Angles)               :  8/ 8 PASS   (1.1s)
------------------------------------------------------------------
  TOTAL                            : 56/56
  Time                             : 617.0s (10.3 min)

  *** ALL PASS ***
==================================================================
```

## Algorithm notes

**BFS orbit enumeration (Theorem B):** Starting from the canonical Hopf
partition, we apply all 120 E₈ Weyl reflections (one per positive root)
to discover new partitions via breadth-first search. Each partition is
encoded as a canonical sorted-tuple-of-sorted-tuples for O(1) hashing.

**Schreier generators (Theorem B):** For each reflection r and each
coset representative gᵢ, compute gᵢ∘r; find the coset representative gⱼ
of the result; the Schreier generator is gⱼ⁻¹∘gᵢ∘r. Closure under
composition yields the full stabilizer.

**Kernel identification (Theorem B):** The 24-element kernel is
identified as 2T = SL(2,F₃) by computing: order distribution
{1:1, 2:1, 3:8, 4:6, 6:8}, center of order 2, derived subgroup of
order 8 that is non-abelian with distribution {1:1, 2:1, 4:6} (= Q₈).
These invariants uniquely determine the binary tetrahedral group.

## License

MIT License — see `LICENSE`.

## Citation

```bibtex
@article{irwin2026kinematic,
  title   = {Kinematic Spinors on Division-Algebra Root Systems:
             A First-Principles Derivation of the Standard Model Weyl Group},
  author  = {Irwin, Klee},
  year    = {2026},
  note    = {Code: \url{https://github.com/QGRKlee/CCT-StandardModel}}
}
```
