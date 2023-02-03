import pandas as pd
import numpy as np

def parse_scores(scores_path):
    scores_data = pd.read_csv(scores_path, sep='\t', header=0)
    nlq = scores_data["-log10q"]
    f10 = np.mean(nlq >= 1)
    f5 = np.mean(nlq >= -np.log10(0.05))
    n = len(nlq)
    
    return scores_data, f10, f5, n

def main(scores_paths, clusters, out_path_data, out_path_summary):
    dfs = []
    summary_records = []
    for sp, c in zip(scores_paths, clusters):
        scores_data, f10, f5, n = parse_scores(sp)

        scores_data['Label'] = c
        dfs.append(scores_data)

        record = {
            "label": c,
            "fdr10": f10,
            "fdr5": f5,
            "npeaks": n,
        }
        summary_records.append(record)

    data_df = pd.concat(dfs)
    summary_df = pd.DataFrame.from_records(summary_records)

    data_df.to_csv(out_path_data, sep="\t", index=False)
    summary_df.to_csv(out_path_summary, sep="\t", index=False)

out_path_data = snakemake.output["data"]
out_path_summary = snakemake.output["summary"]

scores_paths = snakemake.input["scores"]

clusters = snakemake.params["clusters"]

main(scores_paths, clusters, out_path_data, out_path_summary)
    



