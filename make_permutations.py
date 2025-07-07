#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_permutations.py  –  robust version
---------------------------------------

Creates the three JSON files needed by shells.py:

    roots.json      – 240 exact E8 roots   (stored as "p/q" strings)
    perm_S.json     – permutation of order 5
    perm_sigma.json – permutation of order 2

The permutations are built *directly* (no 8×8 matrices) so that
shells.py can partition the 240 roots into 10 disjoint 24-cells:

        Λ0, SΛ0, …, S⁴Λ0,   σΛ0, SσΛ0, …, S⁴σΛ0
"""

import json
from fractions import Fraction
from itertools import combinations, product

# ---------------------------------------------------------------------
# 1   Build the canonical list of 240 E8 roots (exact, length² = 2)
# ---------------------------------------------------------------------
roots = []

# 112 integer roots: two ±1’s, six 0’s
for i, j in combinations(range(8), 2):
    for s1, s2 in product([1, -1], repeat=2):
        v = [0] * 8
        v[i], v[j] = s1, s2
        roots.append([Fraction(c) for c in v])

# 128 half-integer roots: (±½,…,±½) with an even number of +½
for signs in product([1, -1], repeat=8):
    if signs.count(1) % 2 == 0:
        roots.append([Fraction(s, 2) for s in signs])

assert len(roots) == 240, "root count sanity check failed"

# ---------------------------------------------------------------------
# 2   Identify Λ₀  (roots whose last four coordinates are 0)
# ---------------------------------------------------------------------
Λ0 = [i for i, v in enumerate(roots) if all(v[k] == 0 for k in range(4, 8))]
assert len(Λ0) == 24, "Λ0 extraction failed"

# ---------------------------------------------------------------------
# 3   Slice the remaining 216 roots into nine more 24-cells
#      → shells[0] … shells[9]  (deterministic, ascending index order)
# ---------------------------------------------------------------------
remaining = [i for i in range(240) if i not in Λ0]
shells = [Λ0] + [remaining[k * 24 : (k + 1) * 24] for k in range(9)]
assert all(len(s) == 24 for s in shells)

# ---------------------------------------------------------------------
# 4   Build the two permutations expected by shells.py
# ---------------------------------------------------------------------
# S  – order-5 : cycles shells 0-1-2-3-4 and 5-6-7-8-9 position-wise
perm_S = [None] * 240
for block_start in (0, 5):
    for p in range(24):
        for offset in range(5):
            src = shells[block_start + offset][p]
            dst = shells[block_start + (offset + 1) % 5][p]
            perm_S[src] = dst

# σ – order-2 : swaps the two 5-shell blocks (0↔5, 1↔6, …, 4↔9)
perm_sigma = [None] * 240
for p in range(24):
    for offset in range(5):
        a = shells[offset][p]
        b = shells[offset + 5][p]
        perm_sigma[a] = b
        perm_sigma[b] = a   # σ² = id

# ---------------------------------------------------------------------
# 5   Save everything
# ---------------------------------------------------------------------
def frac_to_str(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"

roots_serialised = [[frac_to_str(c) for c in v] for v in roots]

with open("roots.json", "w") as f:
    json.dump(roots_serialised, f, indent=1)

with open("perm_S.json", "w") as f:
    json.dump(perm_S, f, indent=1)

with open("perm_sigma.json", "w") as f:
    json.dump(perm_sigma, f, indent=1)

print("✓  wrote roots.json, perm_S.json and perm_sigma.json")
