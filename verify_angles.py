#!/usr/bin/env python3
"""Theorem E verification: Semi-conformal angle preservation.

Verifies:
  - D4 root system has angle set {60, 90, 120, 180} degrees
  - Projection along a Hopf fiber normal produces two concentric shells
  - Outer (cuboctahedral) shell: 12 points, angles = {60, 90, 120, 180}
  - Inner (octahedral) shell: angles <= {90, 180}
  - No extraneous angles in either shell
  - Result is consistent across all 2-shell projection directions
"""

import sys
import numpy as np
from itertools import combinations
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
print("THEOREM E: Semi-Conformal Angle Preservation")
print("=" * 60)

d4 = build_d4_roots()
N = len(d4)
d4_map = root_index_map(d4)

# ------ D4 angle set in 4D ------
print("\n--- D4 angle set in 4D ---")
d4_angles_4d = set()
for i in range(N):
    for j in range(i + 1, N):
        ni = np.linalg.norm(d4[i])
        nj = np.linalg.norm(d4[j])
        ct = np.dot(d4[i], d4[j]) / (ni * nj)
        ct = np.clip(ct, -1.0, 1.0)
        angle = round(np.degrees(np.arccos(ct)), 1)
        d4_angles_4d.add(angle)

expected_angles = {60.0, 90.0, 120.0, 180.0}
check(d4_angles_4d == expected_angles,
      f"D4 4D angle set = {sorted(d4_angles_4d)} (expected {{60, 90, 120, 180}})")

# ------ Find A2 fibers (partitions of D4 into 4 A2's) ------
print("\n--- Finding A2 fibers ---")
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

# Find one partition
partitions = []
def find_partitions(remaining, current, start):
    if not remaining:
        partitions.append(list(current))
        return
    target = min(remaining)
    for idx in range(start, len(a2_subs)):
        a = a2_subs[idx]
        if target in a and a.issubset(remaining):
            find_partitions(remaining - a, current + [a], idx + 1)

find_partitions(frozenset(range(N)), [], 0)
check(len(partitions) >= 1, f"Found {len(partitions)} D4 partitions into A2's")

# Use the first partition, first fiber
fiber_indices = sorted(list(partitions[0][0]))
fiber_roots = d4[fiber_indices]

# ------ Compute fiber plane and normal space ------
print("\n--- Projection along fiber normal ---")
U, S, Vt = np.linalg.svd(fiber_roots, full_matrices=True)
fiber_plane = Vt[:2]
normal_space = Vt[2:]

def project_to_3d(proj_dir):
    """Project D4 along proj_dir into R^3."""
    n = proj_dir / np.linalg.norm(proj_dir)
    basis_3d = []
    for e in np.eye(4):
        v = e - np.dot(e, n) * n
        for bv in basis_3d:
            v -= np.dot(v, bv) * bv
        norm_v = np.linalg.norm(v)
        if norm_v > 1e-10:
            basis_3d.append(v / norm_v)
        if len(basis_3d) == 3:
            break
    basis_3d = np.array(basis_3d)
    return d4 @ basis_3d.T

def analyze_shells(proj):
    """Analyze shell structure of projected points."""
    norms = np.array([np.linalg.norm(proj[i]) for i in range(N)])
    radii = sorted(set(round(r, 4) for r in norms if r > 1e-10))
    shell_data = {}
    for r in radii:
        shell_idx = [i for i in range(N) if abs(norms[i] - r) < 0.01]
        angles = set()
        for ii in range(len(shell_idx)):
            for jj in range(ii + 1, len(shell_idx)):
                pi_v = proj[shell_idx[ii]]
                pj_v = proj[shell_idx[jj]]
                ni_v = np.linalg.norm(pi_v)
                nj_v = np.linalg.norm(pj_v)
                if ni_v < 1e-10 or nj_v < 1e-10:
                    continue
                ct = np.dot(pi_v, pj_v) / (ni_v * nj_v)
                ct = np.clip(ct, -1.0, 1.0)
                angle = round(np.degrees(np.arccos(ct)), 0)
                angles.add(angle)
        # Count distinct projected points
        pts = set(tuple(np.round(proj[i], 5)) for i in shell_idx)
        shell_data[r] = (len(shell_idx), len(pts), angles)
    return shell_data, radii

# ------ Find 2-shell projection directions ------
# Scan the normal space for directions giving clean 2-shell structure
print("\n--- Scanning for 2-shell projections ---")
d4_angle_set = {60.0, 90.0, 120.0, 180.0}
n_scan = 360
two_shell_dirs = []
two_shell_data = []

for ti in range(n_scan):
    theta = np.pi * ti / n_scan
    proj_dir = np.cos(theta) * normal_space[0] + np.sin(theta) * normal_space[1]
    proj = project_to_3d(proj_dir)
    shell_data, radii = analyze_shells(proj)
    if len(radii) == 2:
        outer_r = radii[-1]
        outer_count, outer_distinct, _ = shell_data[outer_r]
        if outer_distinct == 12:  # 12 DISTINCT projected points = cuboctahedral
            two_shell_dirs.append(theta)
            two_shell_data.append((shell_data, radii))

check(len(two_shell_dirs) >= 1,
      f"Found {len(two_shell_dirs)} clean 2-shell projection directions")

# ------ Detailed verification of first 2-shell projection ------
print("\n--- Cuboctahedral shell analysis ---")
if two_shell_data:
    shell_data_0, radii_0 = two_shell_data[0]

    # Outer (cuboctahedral) shell at largest radius
    outer_r = radii_0[-1]
    outer_count, outer_distinct, outer_angles = shell_data_0[outer_r]

    check(outer_distinct == 12,
          f"Outer shell has {outer_distinct} distinct points (expected 12, cuboctahedral)")

    outer_nonzero = outer_angles - {0.0}
    check(outer_nonzero.issubset(d4_angle_set),
          f"Outer shell angles {sorted(outer_nonzero)} <= {{60,90,120,180}}")
    check(outer_nonzero == d4_angle_set,
          f"Outer shell preserves FULL D4 angle set {{60,90,120,180}}")

    # Inner (octahedral) shell at smaller radius
    inner_r = radii_0[0]
    inner_count, inner_distinct, inner_angles = shell_data_0[inner_r]
    inner_nonzero = inner_angles - {0.0}

    check(inner_nonzero.issubset({90.0, 180.0}),
          f"Inner shell angles {sorted(inner_nonzero)} <= {{90, 180}}")
    check(inner_distinct == 6,
          f"Inner shell has {inner_distinct} distinct points (expected 6, octahedral)")
else:
    for _ in range(5):
        check(False, "Shell analysis", "no 2-shell direction found")

# ------ Verify consistency across all 2-shell directions ------
print("\n--- Cross-checking all 2-shell directions ---")
all_outer_consistent = True
all_inner_consistent = True

for sd, radii in two_shell_data:
    outer_r = radii[-1]
    _, outer_dist, outer_ang = sd[outer_r]
    outer_nz = outer_ang - {0.0}
    if outer_dist != 12 or outer_nz != d4_angle_set:
        all_outer_consistent = False

    inner_r = radii[0]
    _, inner_dist, inner_ang = sd[inner_r]
    inner_nz = inner_ang - {0.0}
    if not inner_nz.issubset({90.0, 180.0}):
        all_inner_consistent = False

check(all_outer_consistent,
      f"Outer shell consistent across all {len(two_shell_dirs)} 2-shell directions")
check(all_inner_consistent,
      f"Inner shell consistent across all {len(two_shell_dirs)} 2-shell directions")

# No extraneous angles in either shell
all_angles_in_any_shell = set()
if two_shell_data:
    sd0 = two_shell_data[0][0]
    for r, (cnt, dist, angles) in sd0.items():
        all_angles_in_any_shell.update(angles - {0.0})
check(all_angles_in_any_shell.issubset(d4_angle_set),
      f"No extraneous angles in any shell (all <= D4 angle set)")

# ================================================================
print("\n" + "=" * 60)
print(f"THEOREM E SUMMARY: {passed}/{passed+failed} claims verified")
if failed > 0:
    print(f"  *** {failed} FAILURES ***")
print("=" * 60)
sys.exit(0 if failed == 0 else 1)
