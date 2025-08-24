from lemmas.lemmas.uniqueness_up_to_conj import prove_uniqueness_up_to_conjugacy
import json

if __name__ == "__main__":
    out = prove_uniqueness_up_to_conjugacy()
    print(json.dumps(out, indent=2, sort_keys=True))
