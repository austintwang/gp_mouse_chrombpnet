import json
import os

def main(folds_dir, out_path):
    in_path = os.paht.join(folds_dir, "fold_0.json")
    with open(in_path) as f:
        data = json.load(f)

    out_data = {"train":[], "valid":[]}
    out_data["test"] = data["test"] + data["train"] + data["valid"]

    with open(out_path, "w") as f:
        json.dump(out_data, f, indent=4)

out_path, = snakemake.output

folds_dir, = snakemake.input

main(folds_dir, out_path)