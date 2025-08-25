# src/runners/a2_hexagon_runner.py
import os, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEM2 = ROOT / "lemmas" / "lemmas"
LEM3 = ROOT / "lemmas" / "lemmas" / "lemmas"
for p in (LEM2, LEM3):
    ps = str(p)
    if ps not in sys.path:
        sys.path.append(ps)

from a2_hexagon_mapping import run

def run_a2_hexagon():
    return run()

if __name__ == "__main__":
    print(json.dumps(run_a2_hexagon(), indent=2, sort_keys=True))
