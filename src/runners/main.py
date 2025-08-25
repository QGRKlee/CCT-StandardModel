# src/runners/main.py
import json

# Import each runner's callable 'run' (all of them already sys.path shim internally)
from runners.geometric_separation_runner import run as run_geo
try:
    from runners.uniqueness_runner import run as run_uniq
except Exception as e:
    run_uniq = None
    _uniq_err = repr(e)

def main():
    results = {}

    # Geometric separation (σ & S commute constraints / shell stats)
    try:
        results["geometric_separation"] = run_geo()
    except Exception as e:
        results["geometric_separation"] = {"ok": False, "error": repr(e)}

    # Uniqueness up to conjugacy (stabilizer BFS certificate)
    if run_uniq:
        try:
            results["uniqueness_up_to_conj"] = run_uniq(max_layers=2, tol=1e-10)
        except Exception as e:
            results["uniqueness_up_to_conj"] = {"ok": False, "error": repr(e)}
    else:
        results["uniqueness_up_to_conj"] = {"ok": False, "error": _uniq_err}

    print(json.dumps(results, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
