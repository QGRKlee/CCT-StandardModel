from lemmas.lemmas.a2_hexagon_mapping import verify_a2_hexagon_mapping
import json

if __name__ == "__main__":
    out = verify_a2_hexagon_mapping()
    print(json.dumps(out, indent=2, sort_keys=True))
