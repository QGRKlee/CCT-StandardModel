# src/runners/main.py
import json

from geometric_separation_runner import run as run_geometric_separation
from uniqueness_runner import run as run_uniqueness

def main():
    results = {
        "geometric_separation": run_geometric_separation(),
        "uniqueness_up_to_conjugacy": run_uniqueness(),
    }
    print(json.dumps(results, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
