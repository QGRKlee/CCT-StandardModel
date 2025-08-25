# src/runners/chirality_runner.py
import os, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEM2 = ROOT / "lemmas" / "lemmas"
LEM3 = ROOT / "lemmas" / "lemmas" / "lemmas"
for p in (LEM2, LEM3):
    ps = str(p)
    if ps not in sys.path:
        sys.path.append(ps)

from chirality_checks import check_chirality

def run():
    """Return chirality diagnostics as a dict."""
    return check_chirality()

if __name__ == "__main__":
    out = run()
    print(json.dumps(out, indent=2, sort_keys=True))
