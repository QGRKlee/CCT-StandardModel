# === [LEMMA-10 A2→SM STUB v1] ===
from __future__ import annotations
import numpy as np

# expected hooks
try:
    from sm_embedding import get_generators
except Exception:
    def get_generators():
        raise NotImplementedError("Provide SM generators (su3, su2, u1Y)")

try:
    from a2_finder import find_A2_subalgebras
except Exception:
    def find_A2_subalgebras():
        raise NotImplementedError("Provide A2 subalgebra finder")

def _comm(A, B):
    return A @ B - B @ A

def _normF(M):
    return float(np.linalg.norm(M, ord="fro"))

def _span_projector(mats: list[np.ndarray]) -> np.ndarray:
    """Build projector onto span{vec(G_i)} using pseudoinverse; works for small 8x8 blocks."""
    if not mats:
        raise ValueError("No matrices provided to build a span projector.")
    G = np.stack([g.reshape(-1) for g in mats], axis=1)  # (64 x m)
    return G @ np.linalg.pinv(G)

def verify_a2_hexagon_mapping() -> dict:
    """
    Verify the A2 hexagon mapping is well‑formed and consistent with labels.
    Replace the body with your actual checks.
    """
    ok = True  # TODO: compute real condition(s)
    return {
        "passed": bool(ok),
        "hexagons_checked": 6,   # TODO
        "labels_ok": True,       # TODO
    }



def check_a2_is_su3_color_compatible(tol: float = 1e-9) -> list[dict]:
    """
    For each A2 found, test closure and overlap with your SU(3)_color span.
    Returns a list of dicts with basic pass/fail diagnostics per A2.
    """
    gens = get_generators()
    su3 = gens["su3"]  # expect list of 8 generators (8x8 each)
    Pspan = _span_projector(su3)

    a2s = find_A2_subalgebras()  # each: {"H1","H2","E_plus","E_minus"}
    out = []
    for idx, a2 in enumerate(a2s):
        H1, H2 = a2["H1"], a2["H2"]
        Ep, Em = a2["E_plus"], a2["E_minus"]

        # su(2)×su(2) style checks inside A2 (Cartan commute; [Ep,Em] ∝ H)
        ok_cartan = _normF(_comm(H1, H2)) <= tol

        # Very light-weight weight checks (placeholders; adapt with your exact normalization):
        # Expect roughly [H1, Ep] ≈ Ep and [H2, Ep] ≈ Ep up to scaling; similarly for Em with minus.
        ok_raise = _normF(_comm(H1, Ep) - Ep) <= 1e-6 and _normF(_comm(H2, Ep) - Ep) <= 1e-6
        ok_lower = _normF(_comm(H1, Em) + Em) <= 1e-6 and _normF(_comm(H2, Em) + Em) <= 1e-6

        # Residual outside su3 span
        def resid(M):
            v = M.reshape(-1, 1)
            r = v - Pspan @ v
            return float(np.linalg.norm(r))
        resid_sum = sum(resid(X) for X in (H1, H2, Ep, Em))

        out.append({
            "a2_index": idx,
            "ok_cartan": bool(ok_cartan),
            "ok_raising_lowering": bool(ok_raise and ok_lower),
            "residual_outside_su3_span": resid_sum,
            "passed": bool(ok_cartan and ok_raise and ok_lower and resid_sum <= 1e-6),
        })
    return out

def check_hypercharge_assignment(tol: float = 1e-9) -> dict:
    """
    Verify that U(1)_Y acts with eigenvalues clustering in the 1:1:1:-3 ratio (up to scale) on your chosen basis.
    Returns the raw spectrum and a coarse ratio diagnostic.
    """
    Y = get_generators()["u1Y"]  # 8x8 matrix
    vals = np.linalg.eigvals(Y)
    vals = np.real_if_close(vals, tol=1e5).astype(float)

    # Sort by absolute value and try to detect 3 nearly-equal vs one ≈ -3× that value
    s = np.sort(vals)
    # crude clustering: average of smallest three vs the remaining extreme
    triplet = np.mean(s[:3])
    extreme = s[-1]
    ratio = (extreme / triplet) if abs(triplet) > tol else float("inf")

    return {
        "eigs_Y_sorted": s.tolist(),
        "triplet_mean": float(triplet),
        "extreme": float(extreme),
        "extreme_over_triplet_mean": float(ratio),
        "approx_1_to_minus3": bool(abs(abs(ratio) - 3.0) <= 0.1),
    }

def run() -> dict:
    """Uniform entry point for runners."""
    return verify_a2_hexagon_mapping()


if __name__ == "__main__":
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))

# === [LEMMA-10 A2→SM STUB v1 END] ===
