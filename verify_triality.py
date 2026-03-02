#!/usr/bin/env python3
"""Theorem D verification: Triality and the 4 = 1 + 3 projection.

Verifies:
  - D4 partitions into 4 disjoint A2 sub-root-systems
  - Exactly 8 such partitions exist
  - The S3 outer automorphism group (triality) acts on these partitions
  - 2 orbits under S3: one of size 2, one of size 6
  - Burnside's lemma check
"""

import sys
import numpy as np
from itertools import combinations, permutations
from e8_utils import build_d4_roots, root_index_map, find_root, TOL

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
print("THEOREM D: Triality and the 4 = 1 + 3 Projection")
print("=" * 60)

# Build D4
d4 = build_d4_roots()
d4_map = root_index_map(d4)
N = len(d4)

check(N == 24, f"D4 has {N} roots (expected 24)")

# ------ Find all A2 sub-root-systems ------
print("\n--- A2 sub-root-systems ---")
a2_subs = []
for combo in combinations(range(N), 6):
    sub = d4[list(combo)]
    if np.linalg.matrix_rank(sub, tol=1e-8) != 2:
        continue
    sub_set = set(combo)
    closed = True
    for i in combo:
        neg_idx = find_root(d4, d4_map, -d4[i])
        if neg_idx not in sub_set:
            closed = False
            break
    if not closed:
        continue
    G = sub @ sub.T
    valid = True
    for ii in range(6):
        for jj in range(ii + 1, 6):
            ip = G[ii, jj]
            if not any(abs(ip - v) < TOL for v in [-2, -1, 1]):
                valid = False
                break
        if not valid:
            break
    if valid:
        a2_subs.append(frozenset(combo))

check(len(a2_subs) == 16, f"Found {len(a2_subs)} A2 sub-root-systems (expected 16)")

# ------ Find all partitions into 4 disjoint A2's ------
print("\n--- Partitions of D4 into 4 disjoint A2's ---")
partitions = []

def find_all_partitions(remaining, current, start):
    if not remaining:
        partitions.append(frozenset(current))
        return
    target = min(remaining)
    for idx in range(start, len(a2_subs)):
        a = a2_subs[idx]
        if target in a and a.issubset(remaining):
            find_all_partitions(remaining - a, current + [a], idx + 1)

find_all_partitions(frozenset(range(N)), [], 0)
check(len(partitions) == 8, f"Found {len(partitions)} partitions (expected 8)")

# ------ Triality matrices ------
print("\n--- S3 triality automorphisms ---")

# Simple roots of D4 (standard labeling):
# a1 = (1,-1,0,0), a2 = (0,1,-1,0), a3 = (0,0,1,-1), a4 = (0,0,1,1)
# Dynkin diagram: a1--a2--a3 with a4 also connected to a2 (central node)
a1 = np.array([1., -1., 0., 0.])
a2 = np.array([0., 1., -1., 0.])
a3 = np.array([0., 0., 1., -1.])
a4 = np.array([0., 0., 1., 1.])
A = np.column_stack([a1, a2, a3, a4])  # (4, 4)
Ai = np.linalg.inv(A)

# tau: outer automorphism of order 3 cycling a1 -> a3 -> a4 -> a1, fixing a2
Bt = np.column_stack([a3, a2, a4, a1])
tau = np.round(Bt @ Ai, 10)

# sigma: outer automorphism of order 2 swapping a1 <-> a3, fixing a2 and a4
# Actually: swapping a3 <-> a4 (the two outer nodes on same side)
Bs = np.column_stack([a1, a2, a4, a3])
sig = np.round(Bs @ Ai, 10)

# Verify orders
tau3 = np.linalg.matrix_power(tau, 3)
check(np.allclose(tau3, np.eye(4), atol=TOL), "tau^3 = I (order 3)")
sig2 = sig @ sig
check(np.allclose(sig2, np.eye(4), atol=TOL), "sigma^2 = I (order 2)")

# tau is NOT a signed permutation matrix (it's an outer automorphism)
is_signed_perm = all(
    sum(abs(tau[i, j]) > 0.1 for j in range(4)) == 1 for i in range(4)
)
check(not is_signed_perm, "tau is NOT a signed permutation (outer automorphism)")

# Verify the explicit matrix from the paper:
# tau = (1/2)[[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[-1,1,1,-1]]
tau_paper = 0.5 * np.array([
    [ 1,  1,  1,  1],
    [ 1,  1, -1, -1],
    [ 1, -1,  1, -1],
    [-1,  1,  1, -1],
])
# tau_paper might be tau or a conjugate; check if it has order 3 and maps D4 to D4
tau_paper_3 = np.linalg.matrix_power(tau_paper, 3)
check(np.allclose(tau_paper_3, np.eye(4), atol=TOL),
      "Paper's tau matrix has order 3")

# Build S3 group
S3_mats = {}
S3_mats["e"] = np.eye(4)
S3_mats["tau"] = tau
S3_mats["tau2"] = tau @ tau
S3_mats["sig"] = sig
S3_mats["sig_tau"] = sig @ tau
S3_mats["sig_tau2"] = sig @ (tau @ tau)

# Verify all 6 elements map D4 roots to D4 roots
all_preserve = True
for name, M in S3_mats.items():
    for i in range(N):
        img = M @ d4[i]
        idx = find_root(d4, d4_map, img)
        if idx < 0:
            all_preserve = False
            break
    if not all_preserve:
        break
check(all_preserve, "All 6 S3 elements map D4 roots to D4 roots")
check(len(S3_mats) == 6, f"S3 has {len(S3_mats)} elements (expected 6)")

# ------ S3 action on partitions ------
print("\n--- S3 action on 8 partitions ---")

def apply_matrix_to_partition(M, partition):
    """Apply a 4x4 matrix to a D4 partition."""
    new_fibers = []
    for fiber in partition:
        new_fiber = set()
        for idx in fiber:
            img = M @ d4[idx]
            new_idx = find_root(d4, d4_map, img)
            assert new_idx >= 0, f"Image of root {idx} under matrix not in D4"
            new_fiber.add(new_idx)
        new_fibers.append(frozenset(new_fiber))
    return frozenset(new_fibers)

# Index partitions
part_list = list(partitions)
part_to_idx = {p: i for i, p in enumerate(part_list)}

# Compute S3 action as permutation of 8 partitions
action_well_defined = True
perm_dict = {}
for name, M in S3_mats.items():
    perm = []
    for p in part_list:
        new_p = apply_matrix_to_partition(M, p)
        if new_p not in part_to_idx:
            action_well_defined = False
            break
        perm.append(part_to_idx[new_p])
    if not action_well_defined:
        break
    perm_dict[name] = tuple(perm)

check(action_well_defined, "S3 action on 8 partitions is well-defined")

# Find orbits
visited = [False] * 8
orbits = []
for i in range(8):
    if visited[i]:
        continue
    orb = set()
    queue = [i]
    while queue:
        x = queue.pop()
        if x in orb:
            continue
        orb.add(x)
        visited[x] = True
        for name, perm in perm_dict.items():
            y = perm[x]
            if y not in orb:
                queue.append(y)
    orbits.append(orb)

orbit_sizes = sorted(len(o) for o in orbits)
check(len(orbits) == 2, f"Number of orbits: {len(orbits)} (expected 2)")
check(orbit_sizes == [2, 6], f"Orbit sizes: {orbit_sizes} (expected [2, 6])")

# Burnside's lemma: sum of |Fix(g)| / |G| = number of orbits
total_fixed = 0
for name, perm in perm_dict.items():
    fixed = sum(1 for i in range(8) if perm[i] == i)
    total_fixed += fixed
burnside = total_fixed / 6
check(abs(burnside - 2.0) < TOL, f"Burnside: sum(|Fix|)/6 = {burnside} (expected 2)")

# ================================================================
print("\n" + "=" * 60)
print(f"THEOREM D SUMMARY: {passed}/{passed+failed} claims verified")
if failed > 0:
    print(f"  *** {failed} FAILURES ***")
print("=" * 60)
sys.exit(0 if failed == 0 else 1)
