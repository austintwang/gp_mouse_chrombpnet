import pandas as pd
import numpy as np

def parse_scores(scores_path):
    scores_data = pd.read_csv(scores_path, sep='\t', header=0)
    nlq = scores_data["prof_nlq"]
    f10 = np.mean(nlq >= 1)
    f1 = np.mean(nlq >= 2)
    n = len(nlq)
    
    return scores_data, f10, f1, n

def parse_counts(counts_path):
    counts_data = pd.read_csv(counts_path, sep='\t', header=0)
    log1pcounts = counts_data["log1p_true_counts"]
    avg_logcounts = np.mean(log1pcounts)

    return counts_data, avg_logcounts

def parse_count_scores(count_scores_path):
    counts_scores = pd.read_csv(count_scores_path, sep='\t', header=0)

    return counts_scores

def main(scores_paths, count_scores_paths, counts_paths, clusters, out_path_data, out_path_summary):
    dfs = []
    summary_records = []
    for sp, csp, cp, c in zip(scores_paths, count_scores_paths, counts_paths, clusters):
        scores_data, f10, f1, n = parse_scores(sp)
        counts_data, avg_logcounts = parse_counts(cp)
        count_scores_data = parse_count_scores(csp)

        data_merged = pd.merge(scores_data, counts_data,  how='left', left_on=['chrom','summit_pos'], right_on=['chrom','summit_pos'])
        data_merged = pd.merge(data_merged, count_scores_data,  how='left', left_on=['chrom','summit_pos'], right_on=['chrom','summit_pos'])
        data_merged['Label'] = c
        dfs.append(data_merged)

        record = {
            "label": c,
            "fdr10": f10,
            "fdr1": f1,
            "npeaks": n,
            "coverage": avg_logcounts
        }
        summary_records.append(record)

    data_df = pd.concat(dfs)
    summary_df = pd.DataFrame.from_records(summary_records)

    data_df.to_csv(out_path_data, sep="\t", index=False)
    summary_df.to_csv(out_path_summary, sep="\t", index=False)

out_path_data = snakemake.output["data"]
out_path_summary = snakemake.output["summary"]

scores_paths = snakemake.input["scores"]
counts_paths = snakemake.input["counts"]
count_scores_paths = snakemake.input["count_scores"]


clusters = snakemake.params["clusters"]

main(scores_paths, count_scores_paths, counts_paths, clusters, out_path_data, out_path_summary)
    



