# src/runners/geometric_separation_runner.py
import os, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEM2 = ROOT / "lemmas" / "lemmas"
LEM3 = ROOT / "lemmas" / "lemmas" / "lemmas"
for p in (LEM2, LEM3):
    ps = str(p)
    if ps not in sys.path:
        sys.path.append(ps)

from geometric_separation import check_geometric_separation

def run():
    """Return the geometric-separation check as a dict."""
    return check_geometric_separation()

if __name__ == "__main__":
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))

