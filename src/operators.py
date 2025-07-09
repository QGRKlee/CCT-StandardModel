import sympy as sp
import json
import os
import numpy as np
from scipy.linalg import expm, logm
from os.path import dirname, abspath, join

# Paths for output
data_dir = abspath(join(dirname(__file__), '..', 'data'))
VECTOR_S_JSON    = join(data_dir, 'S_matrix.json')
VECTOR_SIGMA_JSON= join(data_dir, 'sigma_matrix.json')
SPIN_S_NPY       = join(data_dir, 'S_spin.npy')
SPIN_SIGMA_NPY   = join(data_dir, 'sigma_spin.npy')
os.makedirs(data_dir, exist_ok=True)

# ----- Vector-rep generators -----

def construct_S_vector():
    """
    S: 72° isoclinic rotations (order-5) in planes (0,4),(1,5),(2,6),(3,7)
    """
    theta = 2 * sp.pi / 5
    cos_t, sin_t = sp.cos(theta), sp.sin(theta)
    S = sp.zeros(8, 8)
    for j in range(4):
        i, k = j, j + 4
        S[i, i], S[i, k] = cos_t, -sin_t
        S[k, i], S[k, k] =  sin_t,  cos_t
    return S


def construct_sigma_vector():
    """
    σ: 90° rotations (order-4) in same planes
    """
    phi = sp.pi / 2
    cos_p, sin_p = sp.cos(phi), sp.sin(phi)
    sigma = sp.zeros(8, 8)
    for j in range(4):
        i, k = j, j + 4
        sigma[i, i], sigma[i, k] =  cos_p, -sin_p
        sigma[k, i], sigma[k, k] =  sin_p,  cos_p
    return sigma


def validate_operators():
    """Check orders and commutation of vector cycles"""
    S, sigma = construct_S_vector(), construct_sigma_vector()
    # Order checks
    assert (S**5).equals(sp.eye(8)), "S_vec^5 ≠ I"
    assert (sigma**4).equals(sp.eye(8)), "sigma_vec^4 ≠ I"
    # Commutator
    comm = S * sigma - sigma * S
    assert comm.norm() == 0, f"[S, σ] ≠ 0 (norm = {comm.norm()})"
    print("✅ Vector operators validated: orders and commutation OK")
    return S, sigma


def dump_vector():
    print("Constructing vector representation operators...")
    S, sigma = validate_operators()
    S_np = np.array(S.evalf(), dtype=float)
    sigma_np = np.array(sigma.evalf(), dtype=float)
    json.dump(S_np.tolist(), open(VECTOR_S_JSON, 'w'))
    json.dump(sigma_np.tolist(), open(VECTOR_SIGMA_JSON, 'w'))
    print(f"✅ Wrote {os.path.basename(VECTOR_S_JSON)} and {os.path.basename(VECTOR_SIGMA_JSON)}")

# ----- Spinor-rep generators via half-log mapping -----

def build_spinor():
    """Lift vector-cycle ops to spinor reps with half-angle exponentials"""
    print("Constructing spinor representation operators...")
    S_vec = np.array(json.load(open(VECTOR_S_JSON)), dtype=float)
    sigma_vec = np.array(json.load(open(VECTOR_SIGMA_JSON)), dtype=float)
    L = logm(S_vec)
    Q = logm(sigma_vec)
    S_spin = expm(0.5 * L)
    sigma_spin = expm(0.5 * Q)
    # Check commutation within tolerance
    comm_spin = S_spin @ sigma_spin - sigma_spin @ S_spin
    norm = np.linalg.norm(comm_spin)
    tol = 1e-8
    if norm < tol:
        print(f"✅ Spinor operators commute (norm {norm:.2e} < tol)")
    else:
        raise AssertionError(f"[S_spin, σ_spin] ≠ 0 (norm = {norm})")
    np.save(SPIN_S_NPY, S_spin)
    np.save(SPIN_SIGMA_NPY, sigma_spin)
    print(f"✅ Wrote {os.path.basename(SPIN_S_NPY)} and {os.path.basename(SPIN_SIGMA_NPY)}")

if __name__ == '__main__':
    dump_vector()
    print()
    build_spinor()
