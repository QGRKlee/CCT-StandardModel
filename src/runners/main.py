# src/runners/main.py
# Run available lemma checks and emit a single JSON.

from __future__ import annotations
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]  # .../CCT-StandardModel
SRC  = ROOT / "src"
LEM2 = ROOT / "lemmas" / "lemmas"
LEM3 = ROOT / "lemmas" / "lemmas" / "lemmas"  # your tree shows 3 levels

for p in (SRC, LEM2, LEM3):
    ps = str(p)
    if ps not in sys.path:
        sys.path.append(ps)

# --- imports from src (safe even if not used) ---
try:
    from shells import main as build_shells
except Exception:
    build_shells = None

try:
    from isoclinic import main as verify_isoclinic
except Exception:
    verify_isoclinic = None

# --- lemma imports (use the actual function names in your files) ---
from geometric_separation import check_geometric_separation
from uniqueness_up_to_conj import sample_conjugacy_certificate

def run_all():
    core = {}

    if build_shells:
        try:
            build_shells()
            core["shells"] = "ok"
        except Exception as e:
            core["shells"] = f"error: {e}"

    if verify_isoclinic:
        try:
            verify_isoclinic()
            core["isoclinic"] = "ok"
        except Exception as e:
            core["isoclinic"] = f"error: {e}"

    results = []

    try:
        geo = check_geometric_separation()
        results.append({"lemma": "geometric_separation", **geo})
    except Exception as e:
        results.append({"lemma": "geometric_separation", "error": str(e)})

    try:
        uniq = sample_conjugacy_certificate(max_layers=2)
        results.append({"lemma": "uniqueness_up_to_conj", **uniq})
    except Exception as e:
        results.append({"lemma": "uniqueness_up_to_conj", "error": str(e)})

    payload = {"core": core, "lemmas": results}
    print(json.dumps(payload, indent=2, sort_keys=True))

if __name__ == "__main__":
    run_all()
