# === [LEMMA-8 UNIQUENESS-UP-TO-CONJ STUB v1] ===
from __future__ import annotations
import numpy as np

# expected hooks from repo
try:
    from stabilizer import get_stabilizer_generators_sigma, get_stabilizer_generators_S
except Exception:
    def get_stabilizer_generators_sigma(): return []
    def get_stabilizer_generators_S(): return []

try:
    from isoclinic import get_sigma, get_S
except Exception:
    def get_sigma(): raise NotImplementedError
    def get_S(): raise NotImplementedError

def _charpoly_coeffs(M: np.ndarray) -> np.ndarray:
    vals = np.linalg.eigvals(M)
    return np.sort_complex(vals)

def invariants_pair(S: np.ndarray, sig: np.ndarray) -> dict:
    return {
        "spec_S": _charpoly_coeffs(S),
        "spec_sig": _charpoly_coeffs(sig),
        "spec_Ssig": _charpoly_coeffs(S @ sig),
    }

def _conj(G: np.ndarray, M: np.ndarray) -> np.ndarray:
    return G @ M @ np.linalg.inv(G)

def sample_conjugacy_certificate(num_layers: int = 3, tol: float = 1e-10) -> dict:
    S = get_S()
    sig = get_sigma()
    inv0 = invariants_pair(S, sig)

    Gsig = get_stabilizer_generators_sigma()
    GS   = get_stabilizer_generators_S()
    gens = (Gsig or []) + (GS or [])

    if not gens:
        return {"ok": False, "reason": "No stabilizer generators available"}

    seen = set()
    I = np.eye(S.shape[0])
    frontier = [I]

    for layer in range(1, num_layers + 1):
        new_frontier = []
        for C in frontier:
            for G in gens:
                C2 = G @ C
                h = tuple(np.round(C2.flatten(), 6))
                if h in seen:
                    continue
                seen.add(h)

                Sprime = _conj(C2, S)
                sigprime = _conj(C2, sig)
                inv = invariants_pair(Sprime, sigprime)

                same = all(np.allclose(inv0[k], inv[k], atol=tol) for k in inv0)
                if not same:
                    return {"ok": False, "layer": layer, "bad_invariants": inv}
                new_frontier.append(C2)
        frontier = new_frontier

    return {"ok": True, "tested": len(seen), "layers": num_layers}

if __name__ == "__main__":
    print(sample_conjugacy_certificate())
# === [LEMMA-8 UNIQUENESS-UP-TO-CONJ STUB v1 END] ===
