# === [LEMMA-9 CHIRALITY STUB v1] ===
from __future__ import annotations
import numpy as np

# hooks expected
try:
    from shells import get_shell_orthonormal_bases
except Exception:
    def get_shell_orthonormal_bases(): raise NotImplementedError("Provide orthonormal bases per shell")

try:
    from isoclinic import get_sigma
except Exception:
    def get_sigma(): raise NotImplementedError("Provide σ (order-4)")

def orientation_sign(M: np.ndarray, basis: np.ndarray) -> float:
    """Sign of det of M restricted to span(basis); basis is 8 x d with orthonormal columns."""
    R = basis.T @ (M @ basis)
    det = np.linalg.det(R)
    return float(np.sign(np.real_if_close(det, tol=1e5)))

def check_shell_pairing_chirality(pair_offset: int = 5) -> list[dict]:
    """For each shell k, compare orientation sign on k and k+pair_offset; expect opposite signs."""
    sigma = get_sigma()
    bases = get_shell_orthonormal_bases()
    n = len(bases)
    out = []
    for k in range(n):
        kp = (k + pair_offset) % n
        try:
            s1 = orientation_sign(sigma, bases[k])
            s2 = orientation_sign(sigma, bases[kp])
            ok = (s1 == -s2)
        except np.linalg.LinAlgError:
            s1 = s2 = float("nan")
            ok = False
        out.append({"k": k, "k_pair": kp, "sign_k": s1, "sign_kpair": s2, "opposite": bool(ok)})
    return out
def run() -> dict:
    """Uniform entry point for runners."""
    return check_chirality()

if __name__ == "__main__":
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))

# === [LEMMA-9 CHIRALITY STUB v1 END] ===
