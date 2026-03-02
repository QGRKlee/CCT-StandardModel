#!/usr/bin/env python3
"""Theorem C verification: GUT breaking chain W(D5) ⊃ W(A4) ⊃ W(A2) × W(A1).

Verifies:
  - W(D5) = (Z/2)^4 ⋊ S5 has order 1,920
  - W(A4) = S5 subgroup has order 120
  - [W(D5) : W(A4)] = 16
  - W(A2) × W(A1) = S3 × S2 has order 12
  - [W(A4) : W(A2)×W(A1)] = 10 = C(5,3)
  - All 10 SM embeddings exist
  - Even parity: all W(D5) elements have det = +1 in signed-permutation representation
"""

import sys
import numpy as np
from itertools import permutations, product as iprod, combinations
from e8_utils import compose_perms, inverse_perm, identity_perm, perm_sign

passed = 0
failed = 0

def check(condition, description, detail=""):
    global passed, failed
    if condition:
        print(f"  [PASS] {description}")
        passed += 1
    else:
        print(f"  [FAIL] {description}" + (f" ({detail})" if detail else ""))
        failed += 1

# ================================================================
print("=" * 60)
print("THEOREM C: GUT Breaking Chain")
print("  W(D5) > W(A4) > W(A2) x W(A1)")
print("=" * 60)

# The 10 shells are labeled 0..9, arranged in 5 perpendicular pairs.
# Convention (from w_d5_sm_breaking.py): pairs = {(0,5), (1,6), (2,7), (3,8), (4,9)}
# We label the 5 pairs as P0..P4.
pairs = [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]
pair_idx = {}  # shell -> pair index
for pi, (a, b) in enumerate(pairs):
    pair_idx[a] = pi
    pair_idx[b] = pi

# ------ Build W(D5) = (Z/2)^4 ⋊ S5 ------
print("\n--- W(D5) construction ---")

def make_pair_perm(sigma):
    """Given a permutation sigma of {0,1,2,3,4} (the 5 pairs),
    return the induced permutation of 10 shells."""
    perm = [0] * 10
    for pi in range(5):
        old_a, old_b = pairs[pi]
        new_a, new_b = pairs[sigma[pi]]
        perm[old_a] = new_a
        perm[old_b] = new_b
    return tuple(perm)

def make_flip(flip_vec):
    """Given a binary vector of length 5 indicating which pairs to swap,
    return the permutation of 10 shells."""
    perm = list(range(10))
    for pi in range(5):
        if flip_vec[pi]:
            a, b = pairs[pi]
            perm[a], perm[b] = perm[b], perm[a]
    return tuple(perm)

# S5 part: all permutations of 5 pairs (120 elements)
s5_perms = set()
for sigma in permutations(range(5)):
    s5_perms.add(make_pair_perm(sigma))

# (Z/2)^4 part: even number of within-pair swaps (16 elements)
z2_4_flips = set()
for flip_vec in iprod([0, 1], repeat=5):
    if sum(flip_vec) % 2 == 0:  # Even parity constraint
        z2_4_flips.add(make_flip(flip_vec))

check(len(s5_perms) == 120, f"|S5| = {len(s5_perms)} (expected 120)")
check(len(z2_4_flips) == 16, f"|(Z/2)^4| = {len(z2_4_flips)} (expected 16)")

# Full W(D5): close S5 and (Z/2)^4 under composition
wd5 = set()
for s in s5_perms:
    for f in z2_4_flips:
        wd5.add(compose_perms(f, s))
        wd5.add(compose_perms(s, f))

# Close under composition to be safe
queue = list(wd5)
while queue:
    g = queue.pop()
    for gen in list(s5_perms)[:6] + list(z2_4_flips):  # Use some generators
        for new in [compose_perms(g, gen), compose_perms(gen, g)]:
            if new not in wd5:
                wd5.add(new)
                queue.append(new)

check(len(wd5) == 1920, f"|W(D5)| = {len(wd5)} (expected 1920)")

# ------ W(A4) = S5 subgroup ------
print("\n--- W(A4) = S5 subgroup ---")
wa4 = s5_perms  # Pair permutations only (no within-pair swaps)
check(len(wa4) == 120, f"|W(A4)| = {len(wa4)} (expected 120)")

# Index check
idx_d5_a4 = len(wd5) // len(wa4) if len(wa4) > 0 else 0
check(idx_d5_a4 == 16, f"[W(D5) : W(A4)] = {idx_d5_a4} (expected 16)")

# ------ W(A2) × W(A1) = S3 × S2 ------
print("\n--- W(A2) x W(A1) = S3 x S2 ---")

# For each choice of 3 "color" pairs and 2 "weak" pairs:
# W(A2) = S3 permuting the 3 color pairs
# W(A1) = S2 permuting the 2 weak pairs
# Total order = 6 × 2 = 12

def build_sm_weyl(color_pairs, weak_pairs):
    """Build W(A2) × W(A1) for a given color/weak split."""
    sm = set()
    # S3 on color pairs × S2 on weak pairs
    for cp in permutations(color_pairs):
        for wp in permutations(weak_pairs):
            sigma = [0] * 5
            for i, p in enumerate(cp):
                sigma[p] = color_pairs[i]
            for i, p in enumerate(wp):
                sigma[p] = weak_pairs[i]
            # Actually: sigma[old_pair] = new_pair
            full_sigma = list(range(5))
            for i, c in enumerate(color_pairs):
                full_sigma[c] = cp[i]
            for i, w in enumerate(weak_pairs):
                full_sigma[w] = wp[i]
            sm.add(make_pair_perm(tuple(full_sigma)))
    return sm

# Canonical choice: color = pairs 0,1,2; weak = pairs 3,4
sm_canon = build_sm_weyl([0, 1, 2], [3, 4])
check(len(sm_canon) == 12, f"|W(A2)xW(A1)| = {len(sm_canon)} (expected 12)")

idx_a4_sm = len(wa4) // len(sm_canon) if len(sm_canon) > 0 else 0
check(idx_a4_sm == 10, f"[W(A4) : W(A2)xW(A1)] = {idx_a4_sm} (expected 10)")

# ------ All 10 SM embeddings ------
print("\n--- 10 SM embeddings ---")
all_10_valid = True
count = 0
for color_combo in combinations(range(5), 3):
    weak = [i for i in range(5) if i not in color_combo]
    sm = build_sm_weyl(list(color_combo), weak)
    if len(sm) != 12:
        all_10_valid = False
        break
    # Verify it's a subgroup of W(A4)
    if not sm.issubset(wa4):
        all_10_valid = False
        break
    count += 1

check(count == 10, f"C(5,3) = {count} SM embeddings (expected 10)")
check(all_10_valid, "All 10 SM embeddings are valid order-12 subgroups of W(A4)")

# ------ Even parity ------
print("\n--- Even parity verification ---")

# For W(D5): each element acts as a permutation of 10 shells.
# In the 10-dimensional representation:
#   - Each pair permutation sigma contributes 2 transpositions per crossing -> even parity
#   - Each within-pair swap is 1 transposition, but W(D5) requires an even number -> even parity
# Therefore det(g) = sign(g as 10-perm) = +1 for all g in W(D5).
all_even = True
for g in wd5:
    sign = perm_sign(g)  # sign of g as a permutation of {0,...,9}
    if sign != 1:
        all_even = False
        break

check(all_even, "All 1,920 W(D5) elements have even parity (det = +1)")

# ================================================================
print("\n" + "=" * 60)
print(f"THEOREM C SUMMARY: {passed}/{passed+failed} claims verified")
if failed > 0:
    print(f"  *** {failed} FAILURES ***")
print("=" * 60)
sys.exit(0 if failed == 0 else 1)
