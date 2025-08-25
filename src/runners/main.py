# src/runners/main.py
import json
from geometric_separation_runner import run as geo_run
from uniqueness_runner import run as uniq_run

RUNNERS = {
    "geometric_separation": geo_run,
    "uniqueness": uniq_run,
}

if __name__ == "__main__":
    results = {name: fn() for name, fn in RUNNERS.items()}
    print(json.dumps(results, indent=2, sort_keys=True))
