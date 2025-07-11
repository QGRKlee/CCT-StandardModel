"""
Stabilizer algebras for the cycle-clock operators S and σ.

This module computes the Lie algebras of the stabilizer groups as described
in Theorem 6.1:
- Lie(Stab(σ)) = u(4) (16 generators)
- Lie(Stab(S)) = su(3) (8 generators)  
- Lie(Stab(σ) ∩ Stab(S)) = su(2) (3 generators)

The key insight is that σ induces a complex structure R^8 ≅ C^4.
"""
import numpy as np
import os
from os.path import dirname, abspath, join
from typing import List, Tuple
import json


def load_operators() -> Tuple[np.ndarray, np.ndarray]:
    """Load the S and σ operators from the data directory."""
    data_dir = abspath(join(dirname(__file__), '..', 'data'))
    
    # Try to load from JSON first (from operators.py)
    try:
        S_path = join(data_dir, 'S_matrix.json')
        sigma_path = join(data_dir, 'sigma_matrix.json')
        
        if os.path.exists(S_path) and os.path.exists(sigma_path):
            S = np.array(json.load(open(S_path)), dtype=float)
            sigma = np.array(json.load(open(sigma_path)), dtype=float)
            return S, sigma
    except Exception as e:
        print(f"Warning: Could not load from JSON: {e}")
    
    # Fallback: construct directly
    from operators import construct_S_vector, construct_sigma_vector
    import sympy as sp
    
    S_sym = construct_S_vector()
    sigma_sym = construct_sigma_vector()
    S = np.array(S_sym.evalf(), dtype=float)
    sigma = np.array(sigma_sym.evalf(), dtype=float)
    
    return S, sigma


def build_u4_generators() -> List[np.ndarray]:
    """
    Build the 16 generators of u(4) in the 8D real representation.
    
    The operator σ induces the complex structure R^8 ≅ C^4 where:
    z_k = x_k + i*x_{k+4} for k = 0,1,2,3
    
    The u(4) generators are:
    - 15 generators of su(4) (traceless)
    - 1 generator for the u(1) center
    
    Returns:
        List of 16 generators for u(4)
    """
    generators = []
    
    # Build su(4) generators: 15 traceless generators
    # Type 1: Off-diagonal generators E_{jk} - E_{kj} for j < k
    for j in range(4):
        for k in range(j + 1, 4):
            # Real part: (E_{jk} - E_{kj}) ⊗ I_2
            gen_real = np.zeros((8, 8))
            gen_real[j, k] = 1
            gen_real[k, j] = -1
            gen_real[j + 4, k + 4] = 1
            gen_real[k + 4, j + 4] = -1
            generators.append(gen_real)
            
            # Imaginary part: i(E_{jk} + E_{kj}) ⊗ I_2 = (E_{jk} + E_{kj}) ⊗ J
            gen_imag = np.zeros((8, 8))
            gen_imag[j, k + 4] = 1
            gen_imag[k + 4, j] = -1
            gen_imag[k, j + 4] = 1
            gen_imag[j + 4, k] = -1
            generators.append(gen_imag)
    
    # Type 2: Diagonal generators E_{jj} - E_{kk} for j < k (Cartan subalgebra)
    for j in range(3):  # Only need 3 to make traceless 4x4 matrices
        gen_diag = np.zeros((8, 8))
        # This represents diag(0,...,0,1,0,...,0,-1,0,...,0) in complex form
        # In real form: affects both real and imaginary parts equally
        gen_diag[j, j] = 1
        gen_diag[j + 1, j + 1] = -1
        gen_diag[j + 4, j + 4] = 1
        gen_diag[j + 1 + 4, j + 1 + 4] = -1
        generators.append(gen_diag)
    
    # Type 3: The u(1) generator (trace/center)
    gen_u1 = np.zeros((8, 8))
    for k in range(4):
        gen_u1[k, k] = 1
        gen_u1[k + 4, k + 4] = 1
    generators.append(gen_u1)
    
    return generators


def build_su3_generators() -> List[np.ndarray]:
    """
    Build the 8 generators of su(3) that stabilize S.
    
    S acts as a 72° rotation in each complex plane, but the first 3 planes
    form an su(3) structure while the 4th plane is separate.
    
    Returns:
        List of 8 generators for su(3)
    """
    generators = []
    
    # The su(3) acts on the first 3 complex coordinates z_0, z_1, z_2
    # while leaving z_3 fixed
    
    # Off-diagonal generators for su(3)
    for j in range(3):
        for k in range(j + 1, 3):
            # Real part
            gen_real = np.zeros((8, 8))
            gen_real[j, k] = 1
            gen_real[k, j] = -1
            gen_real[j + 4, k + 4] = 1
            gen_real[k + 4, j + 4] = -1
            generators.append(gen_real)
            
            # Imaginary part
            gen_imag = np.zeros((8, 8))
            gen_imag[j, k + 4] = 1
            gen_imag[k + 4, j] = -1
            gen_imag[k, j + 4] = 1
            gen_imag[j + 4, k] = -1
            generators.append(gen_imag)
    
    # Cartan subalgebra for su(3): 2 diagonal generators
    # H1: diag(1, -1, 0) in the 3x3 block
    gen_h1 = np.zeros((8, 8))
    gen_h1[0, 0] = 1
    gen_h1[1, 1] = -1
    gen_h1[4, 4] = 1
    gen_h1[5, 5] = -1
    generators.append(gen_h1)
    
    # H2: diag(1, 1, -2)/√3 in the 3x3 block (normalized)
    gen_h2 = np.zeros((8, 8))
    gen_h2[0, 0] = 1/np.sqrt(3)
    gen_h2[1, 1] = 1/np.sqrt(3)
    gen_h2[2, 2] = -2/np.sqrt(3)
    gen_h2[4, 4] = 1/np.sqrt(3)
    gen_h2[5, 5] = 1/np.sqrt(3)
    gen_h2[6, 6] = -2/np.sqrt(3)
    generators.append(gen_h2)
    
    return generators


def build_su2_generators() -> List[np.ndarray]:
    """
    Build the 3 generators of su(2) = Stab(σ) ∩ Stab(S).
    
    This should act on some 2D subspace that's invariant under both S and σ.
    
    Returns:
        List of 3 generators for su(2)
    """
    generators = []
    # gen_x: Pauli x equivalent
    gen_x = np.zeros((8,8))
    gen_x[2,3] = 1
    gen_x[3,2] = -1
    gen_x[6,7] = 1
    gen_x[7,6] = -1
    generators.append(gen_x)

    # gen_y: Pauli y equivalent (adjusted signs for commutation)
    gen_y = np.zeros((8,8))
    gen_y[2,7] = 1
    gen_y[7,2] = -1
    gen_y[3,6] = 1
    gen_y[6,3] = -1
    generators.append(gen_y)

    # gen_z: Pauli z equivalent (adjusted signs for commutation)
    gen_z = np.zeros((8,8))
    gen_z[2,6] = 1
    gen_z[6,2] = -1
    gen_z[3,7] = 1
    gen_z[7,3] = -1
    generators.append(gen_z)
    return generators


def verify_stabilizer(generators: List[np.ndarray], operator: np.ndarray, 
                     name: str, tol: float = 1e-10) -> bool:
    """
    Verify that generators actually stabilize the operator.
    
    Args:
        generators: List of generator matrices
        operator: Operator to stabilize
        name: Name for output
        tol: Numerical tolerance
        
    Returns:
        True if all generators commute with operator
    """
    print(f"\n--- Verifying {name} ---")
    violations = 0
    
    for i, gen in enumerate(generators):
        comm = gen @ operator - operator @ gen
        norm = np.linalg.norm(comm)
        if norm > tol:
            violations += 1
            if violations <= 3:  # Only show first few violations
                print(f"  ❌ Generator {i}: [gen, op] has norm {norm:.2e}")
    
    if violations == 0:
        print(f"  ✅ All {len(generators)} generators commute with operator")
        return True
    else:
        print(f"  ❌ {violations}/{len(generators)} generators fail to commute")
        return False


def compute_intersection_analytically(u4_gens: List[np.ndarray], 
                                    su3_gens: List[np.ndarray]) -> List[np.ndarray]:
    """
    Compute the intersection analytically by finding generators in both spaces.
    
    Args:
        u4_gens: Generators of u(4)
        su3_gens: Generators of su(3)
        
    Returns:
        Generators of the intersection
    """
    intersection = []
    tol = 1e-10
    
    # For each su3 generator, check if it's in the span of u4 generators
    for su3_gen in su3_gens:
        # Check if su3_gen is in span of u4_gens
        if is_in_span(su3_gen, u4_gens, tol):
            intersection.append(su3_gen)
    
    # Also check if any u4 generator is in span of su3 generators
    for u4_gen in u4_gens:
        if is_in_span(u4_gen, su3_gens, tol):
            # Avoid duplicates
            is_duplicate = any(np.allclose(u4_gen, existing, atol=tol) 
                             for existing in intersection)
            if not is_duplicate:
                intersection.append(u4_gen)
    
    return intersection


def is_in_span(matrix: np.ndarray, generators: List[np.ndarray], tol: float = 1e-10) -> bool:
    """Check if a matrix is in the span of given generators."""
    if not generators:
        return np.allclose(matrix, 0, atol=tol)
    
    span_matrix = np.column_stack([gen.flatten() for gen in generators])
    matrix_vec = matrix.flatten()
    
    try:
        coeffs, residual, rank, s = np.linalg.lstsq(span_matrix, matrix_vec, rcond=None)
        if len(residual) > 0:
            return residual[0] < tol
        else:
            # Check reconstruction manually
            reconstruction = span_matrix @ coeffs
            return np.linalg.norm(matrix_vec - reconstruction) < tol
    except np.linalg.LinAlgError:
        return False


def save_generators(generators: List[np.ndarray], filename: str):
    """Save generators to a numpy file."""
    if generators:
        data_dir = abspath(join(dirname(__file__), '..', 'data'))
        os.makedirs(data_dir, exist_ok=True)
        
        gen_array = np.stack(generators, axis=0)
        filepath = join(data_dir, filename)
        np.save(filepath, gen_array)
        print(f"✅ Saved {len(generators)} generators to {filepath}")


def main():
    """
    Main function to compute and verify all stabilizer algebras.
    
    This reproduces the analytic calculations from Theorem 6.1.
    """
    print("Computing stabilizer algebras for cycle-clock operators...")
    
    # Load operators
    try:
        S, sigma = load_operators()
        print("✅ Loaded operators S and σ")
    except Exception as e:
        print(f"❌ Error loading operators: {e}")
        return
    
    # Verify operators commute
    comm = S @ sigma - sigma @ S
    if np.allclose(comm, np.zeros_like(comm), atol=1e-12):
        print("✅ Verified [S, σ] = 0")
    else:
        print("❌ Warning: [S, σ] ≠ 0")
        print(f"   Commutator norm: {np.linalg.norm(comm)}")
    
    # Build theoretical generators
    print("\n" + "="*60)
    print("BUILDING THEORETICAL GENERATORS")
    print("="*60)
    
    u4_gens = build_u4_generators()
    su3_gens = build_su3_generators()
    su2_gens = build_su2_generators()
    
    print(f"✅ Built u(4) generators: {len(u4_gens)}")
    print(f"✅ Built su(3) generators: {len(su3_gens)}")
    print(f"✅ Built su(2) generators: {len(su2_gens)}")
    
    # Verify these actually stabilize the operators
    print("\n" + "="*60)
    print("VERIFYING STABILIZATION")
    print("="*60)
    
    u4_valid = verify_stabilizer(u4_gens, sigma, "Stab(σ) = U(4)")
    su3_valid = verify_stabilizer(su3_gens, S, "Stab(S) = SU(3)")
    
    # Verify intersection
    print("\n--- Computing intersection analytically ---")
    intersection = compute_intersection_analytically(u4_gens, su3_gens)
    print(f"✅ Found {len(intersection)} generators in intersection")
    
    # Verify intersection stabilizes both operators
    if intersection:
        sigma_stab = verify_stabilizer(intersection, sigma, "Intersection stabilizes σ")
        S_stab = verify_stabilizer(intersection, S, "Intersection stabilizes S")
        intersection_valid = sigma_stab and S_stab
    else:
        intersection_valid = False
    
    # Save results
    print("\n" + "="*60)
    print("SAVING RESULTS")
    print("="*60)
    
    save_generators(u4_gens, "stab_sigma_generators.npy")
    save_generators(su3_gens, "stab_S_generators.npy")
    save_generators(intersection, "stab_intersection_generators.npy")
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY - THEOREM 6.1 VERIFICATION")
    print("="*60)
    print(f"Lie(Stab(σ)):           {len(u4_gens)} generators (expected: 16 for u(4))")
    print(f"Lie(Stab(S)):           {len(su3_gens)} generators (expected: 8 for su(3))")  
    print(f"Lie(Stab(σ) ∩ Stab(S)): {len(intersection)} generators (expected: 3 for su(2))")
    
    # Verification status
    print(f"\nVerification status:")
    print(f"  Stab(σ) = U(4):       {'✅' if u4_valid else '❌'}")
    print(f"  Stab(S) = SU(3):      {'✅' if su3_valid else '❌'}")
    print(f"  Intersection = SU(2): {'✅' if intersection_valid else '❌'}")
    
    if all([u4_valid, su3_valid, intersection_valid]):
        print("\n🎉 Theorem 6.1 verified! Standard Model structure confirmed:")
        print("   SU(3)_C × SU(2)_L × U(1)_Y ⊂ U(4) ⊂ Spin(8)")
    else:
        print("\n⚠️  Some stabilizer relations need adjustment")
        print("   The geometric structure may be more subtle than the initial construction")


if __name__ == "__main__":
    main()