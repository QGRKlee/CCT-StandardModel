"""
Partition the 240 E8 roots into ten disjoint 24-cells (Λ₀…Λ₉) using the *vector* representation of S and σ.
"""
import os
import numpy as np
import math
from roots import generate_roots


def rotation_matrix(theta: float) -> np.ndarray:
    """Return a 2×2 rotation matrix."""
    return np.array([
        [math.cos(theta), -math.sin(theta)],
        [math.sin(theta),  math.cos(theta)]
    ])


def construct_vector_S() -> np.ndarray:
    """
    72° rotation in the four 2-planes (0,4),(1,5),(2,6),(3,7).
    """
    theta = 2 * math.pi / 5
    R = rotation_matrix(theta)
    M = np.eye(8)
    planes = [(0,4), (1,5), (2,6), (3,7)]
    for i, j in planes:
        M[np.ix_([i,j], [i,j])] = R
    return M


def construct_vector_sigma() -> np.ndarray:
    """
    90° rotation in the four 2-planes (0,4),(1,5),(2,6),(3,7).
    """
    theta = math.pi / 2
    R = rotation_matrix(theta)
    M = np.eye(8)
    planes = [(0,4), (1,5), (2,6), (3,7)]
    for i, j in planes:
        M[np.ix_([i,j], [i,j])] = R
    return M


def partition_shells() -> list[np.ndarray]:
    roots = np.array(generate_roots(), dtype=float)
    S = construct_vector_S()
    sigma = construct_vector_sigma()

    # 1. Λ0: roots with last four coords == 0
    mask0 = np.all(np.isclose(roots[:,4:], 0.0), axis=1)
    shell0 = roots[mask0]
    shells = [shell0]

    # 2. Λ1…Λ4: Sᵏ(Λ0), k=1..4
    for k in range(1, 5):
        shells.append((np.linalg.matrix_power(S, k) @ shell0.T).T)

    # 3. Λ5: σ(Λ0)
    shell5 = (sigma @ shell0.T).T
    shells.append(shell5)

    # 4. Λ6…Λ9: Sᵏ(Λ5), k=1..4
    for k in range(1, 5):
        shells.append((np.linalg.matrix_power(S, k) @ shell5.T).T)

    return shells


def validate_shells(shells: list[np.ndarray]) -> None:
    # Check counts, disjointness, completeness
    if len(shells) != 10:
        raise AssertionError(f"Expected 10 shells, got {len(shells)}")
    # Each shell has 24 vectors
    for i, sh in enumerate(shells):
        if sh.shape[0] != 24:
            raise AssertionError(f"Shell {i} has {sh.shape[0]} vectors, expected 24")
    # Disjoint union
    all_pts = np.vstack(shells)
    if all_pts.shape[0] != 240:
        raise AssertionError(f"Total vectors {all_pts.shape[0]}, expected 240")
    # Unique
    uniq = {tuple(v) for v in map(tuple, all_pts)}
    if len(uniq) != 240:
        raise AssertionError("Shells overlap or missing roots")


def main():
    shells = partition_shells()
    validate_shells(shells)
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
    os.makedirs(data_dir, exist_ok=True)
    # Save each shell to data/ as NPY
    for idx, sh in enumerate(shells):
        np.save(os.path.join(data_dir, f"shell_{idx}.npy"), sh)
    print("Successfully partitioned roots into 10 shells and saved to data/")


if __name__ == "__main__":
    main()

# --- hook used by lemmas: returns 10 orthonormal bases, one per 24-cell shell ---
import numpy as np

def get_shell_orthonormal_bases():
    """
    Return list of orthonormal bases (8 x d_k) for each disjoint 24-cell shell Λ_k.
    If your code builds shells differently, replace the builder call below with your function/variable.
    """
    shells = None
    for fn_name in ("build_shells", "get_shells", "partition_shells", "shells"):
        fn = globals().get(fn_name, None)
        if callable(fn):
            shells = fn()
            break
        if fn_name == "shells" and isinstance(globals().get("shells", None), (list, tuple)):
            shells = globals()["shells"]
            break

    if shells is None:
        raise RuntimeError("No shells found — replace call in get_shell_orthonormal_bases() with your actual builder.")

    bases = []
    for sh in shells:
        A = np.array(sh)
        if A.shape[1] == 8 and A.shape[0] != 8:
            A = A.T
        if A.shape[0] != 8:
            A = A.T
        Q, R = np.linalg.qr(A)
        r = np.linalg.matrix_rank(A)
        bases.append(Q[:, :r])
    return bases
