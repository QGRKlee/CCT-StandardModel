#!/usr/bin/env python3
import json, numpy as np
from fractions import Fraction

# Load the permutations and the string-based roots
perm_S     = json.load(open("perm_S.json"))
perm_sigma = json.load(open("perm_sigma.json"))
roots_str  = json.load(open("roots.json"))

# Parse the fraction strings (e.g., "1/2") into floating-point numbers
roots = [[float(Fraction(c)) for c in v] for v in roots_str]
R = np.array(roots)

# Find the indices of the initial 24-cell (D4 subalgebra)
Λ0 = [i for i,v in enumerate(R) if np.allclose(v[4:],0.0)]
shells=[Λ0]

# Generate 4 more shells by repeatedly applying S
cur=Λ0
for _ in range(1,5):
    cur=[perm_S[i] for i in cur]
    shells.append(cur)

# Generate a new base shell by applying σ to the first shell
Λ5=[perm_sigma[i] for i in Λ0]
shells.append(Λ5)

# Generate the last 4 shells by applying S to the new base shell
cur=Λ5
for _ in range(6,10):
    cur=[perm_S[i] for i in cur]
    shells.append(cur)

# --- Verification ---
flat=[j for Λ in shells for j in Λ]
assert len(shells)==10 and all(len(Λ)==24 for Λ in shells)
assert len(set(flat))==240, f"Error: Found {len(set(flat))} unique roots instead of 240."

print("✓ ten disjoint 24-cells verified.")