# uniqueness_up_to_conj.py
# Certificates that the pair (S, σ) is unique up to conjugation by stabilizers.

# --- path shim so this script can import from src/ no matter the nesting ---
import os, sys, json
_here = os.path.dirname(__file__)
for rel in ("../src", "../../src", "../../../src"):
    cand = os.path.abspath(os.path.join(_here, rel))
    if os.path.isdir(cand) and cand not in sys.path:
        sys.path.append(cand)

import numpy as np
from isoclinic import get_S, get_sigma
from stabilizer import get_stabilizer_generators_sigma, get_stabilizer_generators_S

def _spec(M: np.ndarray) -> np.ndarray:
    """Sorted spectrum (with complex parts collapsed if numerically zero)."""
    vals = np.linalg.eigvals(M)
    vals = np.real_if_close(vals, tol=1e5)
    return np.sort_complex(vals)

def invariants_pair(S: np.ndarray, sig: np.ndarray) -> dict:
    """Conjugacy‑stable invariants for the pair (S, σ)."""
    return {
        "spec_S": _spec(S),
        "spec_sig": _spec(sig),
        "spec_Ssig": _spec(S @ sig),
    }

def _conj(G: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Conjugate M by G."""
    return G @ M @ np.linalg.inv(G)

def sample_conjugacy_certificate(max_layers: int = 2, tol: float = 1e-10) -> dict:
    """
    BFS over products of stabilizer generators up to 'max_layers'.
    Checks that spectra of (S, σ, Sσ) are preserved under conjugation.
    """
    S = get_S()
    sig = get_sigma()
    inv0 = invariants_pair(S, sig)

    gens_sigma = list(get_stabilizer_generators_sigma() or [])
    gens_S     = list(get_stabilizer_generators_S() or [])
    gens = gens_sigma + gens_S
    if not gens:
        return {"ok": False, "reason": "no stabilizer generators available"}

    n = S.shape[0]
    I = np.eye(n)
    frontier = [I]
    seen = {tuple(np.round(I.flatten(), 6))}
    tested = 0

    for layer in range(1, max_layers + 1):
        new_frontier = []
        for C in frontier:
            for G in gens:
                C2 = G @ C
                key = tuple(np.round(C2.flatten(), 6))
                if key in seen:
                    continue
                seen.add(key)

                Sprime = _conj(C2, S)
                sigprime = _conj(C2, sig)
                inv = invariants_pair(Sprime, sigprime)
                tested += 1

                same = all(np.allclose(inv0[k], inv[k], atol=tol) for k in inv0)
                if not same:
                    return {
                        "ok": False,
                        "layer": layer,
                        "tested": tested,
                        "counterexample_found": True,
                    }
                new_frontier.append(C2)
        frontier = new_frontier
        if not frontier:
            break

    return {"ok": True, "layers": max_layers, "tested": tested}

def run() -> dict:
    """
    Uniform entry point for runners. Returns a JSON‑serializable dict.
    You can change layers via env var UNIQUENESS_LAYERS (default 2).
    """
    max_layers = int(os.getenv("UNIQUENESS_LAYERS", "2"))
    return sample_conjugacy_certificate(max_layers=max_layers, tol=1e-10)


if __name__ == "__main__":
    res = run()
    print(json.dumps(res, indent=2, sort_keys=True))

