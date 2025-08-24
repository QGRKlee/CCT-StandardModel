# runners/main.py

import json
from pathlib import Path

# Always works (you already have this file)
from .geometric_separation_runner import run as run_geometric_separation

# Add these later, after we create the files
# from .uniqueness_runner import run as run_uniqueness
# from .a2_hexagon_runner import run as run_a2_hexagon

def main():
    out = {}

    # 1) geometric separation
    out["geometric_separation"] = run_geometric_separation()

    # 2) uniqueness up to conjugacy (uncomment after file exists)
    # out["uniqueness_up_to_conj"] = run_uniqueness()

    # 3) A2 hexagon mapping (uncomment after file exists)
    # out["a2_hexagon_mapping"] = run_a2_hexagon()

    # write a single JSON blob next to this file (../data/lemmas_results.json)
    data_dir = Path(__file__).resolve().parent.parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "lemmas_results.json").write_text(json.dumps(out, indent=2))
    print("✅ wrote", data_dir / "lemmas_results.json")

if __name__ == "__main__":
    main()
