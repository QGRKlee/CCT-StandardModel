from lemmas.lemmas.geometric_separation import check_geometric_separation
import json

if __name__ == "__main__":
    out = check_geometric_separation()
    print(json.dumps(out, indent=2, sort_keys=True))
