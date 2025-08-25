# src/runners/main.py
import json
from importlib import import_module

# Map task keys to (module, function)
TASKS = {
    "geometric_separation": ("src.runners.geometric_separation_runner", "run"),
    "uniqueness": ("src.runners.uniqueness_runner", "run"),
    "chirality": ("src.runners.chirality_runner", "run"),
}

def run_all():
    out = {}
    for key, (mod, fn) in TASKS.items():
        m = import_module(mod)
        res = getattr(m, fn)()
        # normalize typical fields
        if isinstance(res, dict):
            res.setdefault("ok", True if "passed" in res and res["passed"] else res.get("ok", None))
        out[key] = res
    return out

if __name__ == "__main__":
    print(json.dumps(run_all(), indent=2, sort_keys=True))
