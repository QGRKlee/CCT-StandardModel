from lemmas.lemmas.chirality_checks import run_chirality_checks
import json

if __name__ == "__main__":
    out = run_chirality_checks()
    print(json.dumps(out, indent=2, sort_keys=True))
