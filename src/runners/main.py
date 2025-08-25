# src/runners/main.py
import json

from geometric_separation_runner import run as run_geometric_separation
from uniqueness_runner import run as run_uniqueness
from chirality_runner import run_chirality
from a2_hexagon_runner import run_a2_hexagon

def main():
    results = {
        "geometric_separation": run_geometric_separation(),
        "uniqueness_up_to_conjugacy": run_uniqueness(),
        "chirality": run_chirality(),
        "a2_hexagon_mapping": run_a2_hexagon(),
    }
    print(json.dumps(results, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
