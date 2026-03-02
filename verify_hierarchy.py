#!/usr/bin/env python3
"""Theorem A verification: Compositional hierarchy A1 ⊂ A2 ⊂ D4 ⊂ E8.

Verifies:
  - Root system construction at each level
  - Kinematic spinor (Level 1): Spin(2) rotor generates A2 from A1
  - Spinorial signature: spinor orbit has 2x elements vs SO image
  - Compositional structure: each level is a compound of the previous
  - Hopf fibration structure at Levels 2 and 3
"""

import sys
import numpy as np
from itertools import combinations
from e8_utils import (build_e8_roots, build_d4_roots, build_a2_roots, build_a1_roots,
                      root_index_map, find_root, verify_root_system,
                      cluster_by_hopf, verify_d4_cartan, TOL)

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
print("THEOREM A: Compositional Hierarchy A1 < A2 < D4 < E8")
print("=" * 60)

# ------ Level 0: A1 ------
print("\n--- Level 0: A1 ---")
a1 = build_a1_roots()
check(a1.shape == (2, 1), "A1 has 2 roots in R^1")
check(np.allclose(np.sum(a1**2, axis=1), 2.0), "A1 roots have norm^2 = 2")

# ------ Level 1: A2 via kinematic spinor ------
print("\n--- Level 1: A2 via Spin(2) rotor ---")
a2 = build_a2_roots()
ok, msg = verify_root_system(a2, expected_count=6, expected_rank=2)
check(ok, "A2 is a valid root system (6 roots, rank 2)", msg)

# The Spin(2) rotor R = exp(pi*e1e2/3) acts as rotation by pi/3 in R^2
theta = np.pi / 3.0
R = np.array([[np.cos(theta), -np.sin(theta)],
              [np.sin(theta),  np.cos(theta)]])

# Check order
Rk = np.eye(2)
order = 0
for k in range(1, 13):
    Rk = Rk @ R
    if np.allclose(Rk, np.eye(2), atol=TOL):
        order = k
        break
check(order == 6, f"Rotor R has order {order} (expected 6)")

# Spinorial signature: R^3 = -I
R3 = np.linalg.matrix_power(R, 3)
check(np.allclose(R3, -np.eye(2), atol=TOL), "R^3 = -I (spinorial signature)")

# Orbit of a seed vector under R generates A2
seed = a2[0]  # Take first A2 root as seed
orbit_spinor = set()
v = seed.copy()
a2_map = root_index_map(a2)
for k in range(6):
    idx = find_root(a2, a2_map, v)
    orbit_spinor.add(idx)
    v = R @ v
check(len(orbit_spinor) == 6, f"Spinor orbit has {len(orbit_spinor)} elements (expected 6)")
check(orbit_spinor == set(range(6)), "Orbit generates all 6 A2 roots")

# SO image: R^2 acts as rotation by 2*pi/3, giving orbit of 3
R_so = np.linalg.matrix_power(R, 2)  # SO(2) image has half the order
orbit_so = set()
v = seed.copy()
for k in range(3):
    idx = find_root(a2, a2_map, v)
    orbit_so.add(idx)
    v = R_so @ v
check(len(orbit_so) == 3, f"SO orbit has {len(orbit_so)} elements (expected 3)")
check(len(orbit_spinor) == 2 * len(orbit_so), "Spinor/SO orbit ratio = 2 (double cover)")

# A2 = 3 x A1 (three antipodal pairs)
pair_count = 0
used = set()
for i in range(6):
    if i in used:
        continue
    for j in range(i+1, 6):
        if j in used:
            continue
        if np.allclose(a2[i], -a2[j], atol=TOL):
            pair_count += 1
            used.update([i, j])
            break
check(pair_count == 3, f"A2 = 3 antipodal pairs (3 x A1)")

# ------ Level 2: D4 ------
print("\n--- Level 2: D4 ---")
d4 = build_d4_roots()
ok, msg = verify_root_system(d4, expected_count=24, expected_rank=4)
check(ok, "D4 is a valid root system (24 roots, rank 4)", msg)

# D4 = 4 x A2: find at least one partition into 4 disjoint A2 subsystems
d4_map = root_index_map(d4)
# Find all A2 subsystems (6-root subsets that are valid A2)
a2_subs = []
for combo in combinations(range(24), 6):
    sub = d4[list(combo)]
    # Quick checks: rank 2, closure under negation, inner products in {-2,-1,1}
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
        for jj in range(ii+1, 6):
            ip = G[ii, jj]
            if not any(abs(ip - v) < TOL for v in [-2, -1, 1]):
                valid = False
                break
        if not valid:
            break
    if valid:
        a2_subs.append(frozenset(combo))

# Find partition: 4 disjoint A2's covering all 24 roots
partition_found = False
def find_partition(remaining, current, a2_list, start):
    global partition_found
    if partition_found:
        return
    if not remaining:
        partition_found = True
        return
    target = min(remaining)
    for idx in range(start, len(a2_list)):
        a = a2_list[idx]
        if target in a and a.issubset(remaining):
            find_partition(remaining - a, current + [a], a2_list, idx + 1)
            if partition_found:
                return

find_partition(frozenset(range(24)), [], list(a2_subs), 0)
check(partition_found, "D4 = 4 x A2 (partition into 4 disjoint A2 subsystems exists)")

# ------ Level 3: E8 ------
print("\n--- Level 3: E8 ---")
e8 = build_e8_roots()
# Count types
integer_count = sum(1 for r in e8 if all(abs(x - round(x)) < TOL for x in r))
half_int_count = len(e8) - integer_count
check(len(e8) == 240, f"E8 has {len(e8)} roots (expected 240)")
check(integer_count == 112, f"Integer-type roots: {integer_count} (expected 112)")
check(half_int_count == 128, f"Half-integer-type roots: {half_int_count} (expected 128)")

ok, msg = verify_root_system(e8, expected_count=240, expected_rank=8)
check(ok, "E8 is a valid root system (240 roots, rank 8)", msg)

# E8 = 10 x D4 via quaternionic Hopf
shells = cluster_by_hopf(e8)
check(len(shells) == 10, f"Hopf partition: {len(shells)} clusters (expected 10)")
all_d4 = True
for i, sh in enumerate(shells):
    sub = e8[sh]
    ok_sub, msg_sub = verify_root_system(sub, expected_count=24, expected_rank=4)
    if not ok_sub:
        all_d4 = False
        break
    ok_cartan, msg_cartan = verify_d4_cartan(sub)
    if not ok_cartan:
        all_d4 = False
        break
check(all_d4, "All 10 Hopf clusters are valid D4 root systems (Cartan matrix verified)")

# ================================================================
print("\n" + "=" * 60)
print(f"THEOREM A SUMMARY: {passed}/{passed+failed} claims verified")
if failed > 0:
    print(f"  *** {failed} FAILURES ***")
print("=" * 60)
sys.exit(0 if failed == 0 else 1)
