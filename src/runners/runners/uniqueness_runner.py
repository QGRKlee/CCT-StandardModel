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

from uniqueness_up_to_conj import sample_conjugacy_certificate

if __name__ == "__main__":
    out = sample_conjugacy_certificate(max_layers=2)
    print(json.dumps(out, indent=2, sort_keys=True))
