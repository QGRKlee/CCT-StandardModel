# Cyclic Group Actions on Division‑Algebra Lattices (Spin(8) → Spin(3))

**Purpose.** This repository collects *machine‑verifiable* certificates and code for the paper:

**“Cyclic Group Actions on Division‑Algebra Lattices: Deriving the Standard Model from Spin(8) → Spin(3)”**

It organizes the proofs, data, and viewers used to certify theorems/lemmas about:
- nested root systems: A1 ⊂ A2 ⊂ D4 ⊂ E8,
- order‑2/3/4/5 clock operators on E8 (pointer/clicker),
- stabilizer algebras inside Spin(8) (U(4), su(3), su(2), SU(2)^4),
- semi‑conformal 8→3 projection invariants (Fibonacci Icosagrid / FIG),
- angle spectra preservation and bipartite structures (4+4, 6+6, 5+5),
- Spin(8) → Spin(3) spinor transport (Pauli/Dirac modules).

The repo contains: (i) **certificates** (human‑readable + machine logs), (ii) **code** (Python; HTML viewer), and (iii) a **Word‑safe draft** of the paper (ASCII symbols only).


## Repo layout

```
.
├── README.md
├── LICENSE
├── .gitignore
├── paper/
│   └── draft_word_safe.md
├── proofs/
│   ├── index.md
│   ├── CERT-01_A1-A2-D4-E8_units/
│   │   └── README.md
│   ├── CERT-02_E8_into_ten_D4_shells/
│   │   └── README.md
│   ├── CERT-03_order5_pointer_two_5cycles/
│   │   └── README.md
│   ├── CERT-04_order4_clicker_Stab_is_U4_dim16/
│   │   └── README.md
│   ├── CERT-05_order2_swap_Stab_is_SU2x4/
│   │   └── README.md
│   ├── CERT-06_intersection_Stab_s_sigma_is_su2_dim3/
│   │   └── README.md
│   ├── CERT-07_angle_spectra_E8_D4_Eisenstein/
│   │   └── README.md
│   ├── CERT-08_FIG_4G_angle_spectrum_match/
│   │   └── README.md
│   ├── CERT-09_semi_conformal_projection_equivariance/
│   │   └── README.md
│   ├── CERT-10_Pauli_to_Dirac_pairing_in_Spin3/
│   │   └── README.md
│   ├── CERT-11_order2_Z_to_Eisenstein_subclock/
│   │   └── README.md
│   ├── CERT-12_SM_chain_SU3xSU2xU1_in_U4_in_Spin8/
│   │   └── README.md
│   ├── TODO-REM-01_order3_operator_u_su3_centralizer/
│   │   └── SPEC.md
│   ├── TODO-REM-02_simultaneous_stabilizers_u_s_sigma/
│   │   └── SPEC.md
│   ├── TODO-REM-03_equivariant_projection_certificate/
│   │   └── SPEC.md
│   ├── TODO-REM-04_russian_doll_hulls_from_true_FIG/
│   │   └── SPEC.md
│   └── TODO-REM-05_hypercharge_quantization_via_C5/
│       └── SPEC.md
└── code/
    ├── requirements.txt
    ├── python/
    │   ├── e8_roots.py
    │   ├── order_operators.py
    │   ├── stabilizers.py
    │   ├── angle_spectra.py
    │   └── __init__.py
    └── viewers/
        └── nested_viewer_fixed.html
```

- **paper/** – ASCII/Word‑safe draft (we will later switch to LaTeX).
- **proofs/** – one subfolder per certificate. `CERT-*` folders hold *ready* proofs; `TODO-REM-*` are the **five** remaining certifications to code.
- **code/** – Python utilities + an interactive HTML viewer.

## Reproducibility quickstart

```bash
# Python 3.10+ recommended
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r code/requirements.txt
# (Planned) run specific proof scripts, e.g.:
# python -m code.python.stabilizers --proof CERT-04
```

### HTML viewer
Open **code/viewers/nested_viewer_fixed.html** in any modern browser (double click). It renders the cuboctahedron hull with two concentric octahedra and runs strict half‑space checks.

## How to publish to GitHub

**Option A — Web UI (fastest)**
1. Create a new GitHub repo (e.g., `cct-sm-spin8`). Leave it empty.
2. Click **“Add file → Upload files”** and drag the *contents* of this folder in.
3. Commit to `main` (or `master`).

**Option B — Command line**
```bash
cd /path/to/cct-sm-spin8-repo
git init
git add .
git commit -m "Initial commit: proofs scaffolding + viewer"
git branch -M main
git remote add origin https://github.com/YOUR-ORG-OR-USER/cct-sm-spin8.git
git push -u origin main
```

After you push, **please send me the GitHub link**. I’ll wire it into the draft and the “Data/Code Availability” section.

## License
See `LICENSE` (placeholder). We can switch to Apache‑2.0/MIT on your instruction.

