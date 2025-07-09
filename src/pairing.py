"""
Compute the unique perpendicular partner for each 24-cell shell.
"""
import os
import numpy as np


def load_shell(i: int) -> np.ndarray:
    """Load shell_i from data/"""
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", f"shell_{i}.npy"))
    return np.load(path)


def find_perpendicular_partners() -> dict[int, int]:
    """
    For each shell Λ_i, find the unique Λ_j (j ≠ i) such that every vector in Λ_i is orthogonal to every vector in Λ_j.
    Returns a dict mapping i→j.
    """
    shells = [load_shell(i) for i in range(10)]
    partners: dict[int, int] = {}

    # tolerance for numerical orthogonality
    tol = 1e-3
    for i, sh_i in enumerate(shells):
        found = False
        for j, sh_j in enumerate(shells):
            if i == j:
                continue
            # compute all pairwise dot-products
            dots = sh_i @ sh_j.T  # (24,24)
            max_dot = np.max(np.abs(dots))
            # treat as perpendicular if all dot products are below tol
            if max_dot < tol:
                partners[i] = j
                found = True
                break
        if not found:
            raise AssertionError(f"No perpendicular partner found for shell {i} (max dot {max_dot})")
    return partners


def main():
    partners = find_perpendicular_partners()
    for i, j in sorted(partners.items()):
        print(f"Shell {i} ⟂ Shell {j}")
    # Save mapping
    out_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "partners.npy"))
    np.save(out_path, partners)
    print(f"Successfully computed and saved shell pairing to {out_path}")


if __name__ == "__main__":
    main()