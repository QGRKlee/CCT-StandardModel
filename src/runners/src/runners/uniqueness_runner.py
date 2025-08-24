# src/runners/uniqueness_runner.py
import os, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEM2 = ROOT / "lemmas" / "lemmas"
LEM3 = ROOT / "lemmas" / "lemmas" / "lemmas"
for p in (LEM2, LEM3):
    ps = str(p)
    if ps not in sys.path:
        sys.path.append(ps)

from uniqueness_up_to_conj import check_uniqueness_up_to_conjugacy

def run():
    """Return the uniqueness-up-to-conjugacy check as a dict."""
    return check_uniqueness_up_to_conjugacy()

if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
