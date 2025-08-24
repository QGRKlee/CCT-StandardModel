# runners/main.py
# Tiny runner that executes our lemma checks and writes a single JSON report.

from __future__ import annotations
import json, os, sys
from pathlib import Path

# --- repo roots / import paths
ROOT = Path(__file__).resolve().parents[2]  # …/CCT-StandardModel
SRC  = ROOT / "src"
LEM  = ROOT / "lemmas" / "lemmas"
DATA = ROOT / "data"
for p in (str(SRC), str(LEM)):
    if p not in sys.path:
        sys.path.append(p)

# --- imports from our codebase
# core pieces (safe to import even if not used in every run)
import numpy as np  # noqa: F401
from shells import main as build_shells          # src/shells.py (generates shell_*.npy)
from isoclinic import main as verify_isoclinic   # src/isoclinic.py (saves R_60_isoclinic.npy)

# lemma runners (these are the stubs we added)
from geometric_separation import check_geometric_separation
from uniqueness_up_to_conj import run_uniqueness_search
from chirality_checks import run_chirality_suite
from a2_hexagon_mapping import run_a2_hexagon_map

def ensure_data():
    DATA.mkdir(exist_ok=True)

def run_all():
    ensure_data()

    # 0) make sure core artifacts exist
    core = {}
    try:
        build_shells()          # writes shell_0.npy … shell_9.npy
        core["shells"] = "ok"
    except Exception as e:
        core["shells"] = f"error: {e}"

    try:
        verify_isoclinic()      # writes R_60_isoclinic.npy (or reuses perpendicular_pairs)
        core["isoclinic"] = "ok"
    except Exception as e:
        core["isoclinic"] = f"error: {e}"

    # 1) lemmas
    results = []

    try:
        geo = check_geometric_separation()
        results.append({"lemma": "geometric_separation", **geo})
    except Exception as e:
        results.append({"lemma": "geometric_separation", "error": str(e)})

    try:
        uniq = run_uniqueness_search(limit=2000)  # small cap so it finishes quickly
        results.append({"lemma": "uniqueness_up_to_conj", **uniq})
    except Exception as e:
        results.append({"lemma": "uniqueness_up_to_conj", "error": str(e)})

    try:
        chir = run_chirality_suite()
        results.append({"lemma": "chirality_checks", **chir})
    except Exception as e:
        results.append({"lemma": "chirality_checks", "error": str(e)})

    try:
        a2 = run_a2_hexagon_map()
        results.append({"lemma": "a2_hexagon_mapping", **a2})
    except Exception as e:
        results.append({"lemma": "a2_hexagon_mapping", "error": str(e)})

    # 2) persist + print
    payload = {"core": core, "lemmas": results}
    out_path = DATA / "lemma_results.json"
    with out_path.open("w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    print(json.dumps(payload, indent=2, sort_keys=True))

if __name__ == "__main__":
    run_all()
