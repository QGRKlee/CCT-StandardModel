"""
Generate and validate the 240 E₈ roots, then write them to data/roots.json.
Corrected to match the paper's specifications.
"""
from fractions import Fraction
from itertools import product, combinations
import json
import os


def generate_type_I():
    """
    Type I roots: two entries ±1 and six zeros, length-squared = 2.
    These are permutations of (±1, 0, 0, 0, 0, 0, 0, 0) with exactly two non-zero entries.
    Total count: C(8,2) * 2^2 = 28 * 4 = 112.
    """
    roots = []
    # Select all unordered pairs of distinct indices for the two non-zero positions
    for i, j in combinations(range(8), 2):
        # All four sign combinations for the two non-zero entries
        for s1 in (Fraction(1), Fraction(-1)):
            for s2 in (Fraction(1), Fraction(-1)):
                vec = [Fraction(0)] * 8
                vec[i] = s1
                vec[j] = s2
                roots.append(tuple(vec))
    return roots


def generate_type_II():
    """
    Type II roots: (1/2)(±1, ±1, ±1, ±1, ±1, ±1, ±1, ±1) with even number of minus signs.
    Total count: 2^7 = 128 (since we fix the constraint of even number of minus signs).
    """
    roots = []
    half = Fraction(1, 2)
    
    # Generate all combinations with even number of minus signs
    for signs in product((half, -half), repeat=8):
        # Count negative signs
        neg_count = sum(1 for s in signs if s < 0)
        if neg_count % 2 == 0:  # Even number of minus signs
            roots.append(tuple(signs))
    
    return roots


def generate_roots():
    """Generate all 240 E₈ roots."""
    type_I = generate_type_I()
    type_II = generate_type_II()
    
    print(f"Type I roots: {len(type_I)}")
    print(f"Type II roots: {len(type_II)}")
    
    roots = type_I + type_II
    
    if len(roots) != 240:
        raise ValueError(f"Expected 240 roots, got {len(roots)}")
    
    return roots


def validate_roots(roots):
    """
    Validate that all roots have the correct properties for E₈:
    - All roots have length √2 (norm = 2)
    - Inner products between distinct roots are in {-2, -1, 0, 1, 2}
    """
    print("Validating root system properties...")
    
    # Check norms
    for i, root in enumerate(roots):
        norm_squared = sum(x * x for x in root)
        if norm_squared != 2:
            raise AssertionError(f"Root {i} has norm² = {norm_squared}, expected 2")
    
    # Check inner products
    allowed_inner_products = {-2, -1, 0, 1, 2}
    inner_product_counts = {ip: 0 for ip in allowed_inner_products}
    
    for i, u in enumerate(roots):
        for j, v in enumerate(roots):
            if i != j:  # Don't check self inner product
                inner_product = sum(x * y for x, y in zip(u, v))
                if inner_product not in allowed_inner_products:
                    raise AssertionError(f"Invalid inner product {inner_product} between roots {i} and {j}")
                inner_product_counts[inner_product] += 1
    
    print("Inner product distribution:")
    for ip, count in inner_product_counts.items():
        print(f"  {ip}: {count} pairs")
    
    print("Root system validation passed!")


def verify_d4_subsystems():
    """
    Verify that choosing any four coordinate axes yields a D₄ subsystem with 24 roots.
    This supports the paper's claim about 24-cells.
    """
    print("\nVerifying D₄ subsystem property...")
    
    roots = generate_roots()
    
    # Test a few different 4-coordinate selections
    test_coords = [
        (0, 1, 2, 3),  # First four coordinates
        (4, 5, 6, 7),  # Last four coordinates
        (0, 2, 4, 6),  # Even coordinates
        (1, 3, 5, 7),  # Odd coordinates
    ]
    
    for coords in test_coords:
        # Count roots that are zero outside the selected coordinates
        d4_roots = []
        for root in roots:
            # Check if root is zero outside selected coordinates
            is_d4_root = all(root[i] == 0 for i in range(8) if i not in coords)
            if is_d4_root:
                # Extract the 4D projection
                projected = tuple(root[i] for i in coords)
                d4_roots.append(projected)
        
        print(f"Coordinates {coords}: {len(d4_roots)} roots form D₄ subsystem")
        
        # Should have exactly 24 roots for D₄
        if len(d4_roots) != 24:
            print(f"  Warning: Expected 24 roots for D₄, got {len(d4_roots)}")
        else:
            print(f"  ✓ Correct D₄ subsystem size")


def main():
    """Main function to generate, validate, and save E₈ roots."""
    print("Generating E₈ root system...")
    
    roots = generate_roots()
    validate_roots(roots)
    verify_d4_subsystems()
    
    # Prepare JSON-serializable list
    json_roots = [[str(c) for c in vec] for vec in roots]
    
    # Ensure data directory exists
    os.makedirs("data", exist_ok=True)
    out_path = os.path.join("data", "roots.json")
    
    # Save to file
    with open(out_path, "w") as f:
        json.dump(json_roots, f, indent=2)
    
    print(f"\nSuccessfully wrote {len(roots)} E₈ roots to {out_path}")
    print("Root generation complete and validated!")


if __name__ == "__main__":
    main()