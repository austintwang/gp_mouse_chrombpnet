import os
import gzip
import numpy as np
import h5py
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import fdrcorrection

def get_powerlaw_at_distance(distances, gamma, min_distance=5000, scale=None):
    """
    From https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/95bde8fff18979e69f7477b2fe240cc4567ff36f/src/hic.py#L185
    """

    assert(gamma > 0)

    #The powerlaw is computed for distances > 5kb. We don't know what the contact freq looks like at < 5kb.
    #So just assume that everything at < 5kb is equal to 5kb.
    #TO DO: get more accurate powerlaw at < 5kb
    distances = np.clip(distances, min_distance, np.Inf)
    log_dists = np.log(distances + 1)

    #Determine scale parameter
    #A powerlaw distribution has two parameters: the exponent and the minimum domain value 
    #In our case, the minimum domain value is always constant (equal to 1 HiC bin) so there should only be 1 parameter
    #The current fitting approach does a linear regression in log-log space which produces both a slope (gamma) and a intercept (scale)
    #Empirically there is a linear relationship between these parameters (which makes sense since we expect only a single parameter distribution)
    #It should be possible to analytically solve for scale using gamma. But this doesn't quite work since the hic data does not actually follow a power-law
    #So could pass in the scale parameter explicity here. Or just kludge it as I'm doing now
    #TO DO: Eventually the pseudocount should be replaced with a more appropriate smoothing procedure.

    #4.80 and 11.63 come from a linear regression of scale on gamma across 20 hic cell types at 5kb resolution. Do the params change across resolutions?
    if scale is None:
        scale = -4.80 + 11.63 * gamma

    powerlaw_contact = np.exp(scale + -1*gamma * log_dists)

    return(powerlaw_contact)

def load_logcts(hdf5_paths_ss, hdf5_paths_xs, true_counts_path):
    chroms = None
    coords = None
    data_ss = []
    for p in hdf5_paths_ss:    
        with h5py.File(p) as f:
            indicators = f["coords/coords_peak"][:]
            logcts = f["predictions/logcounts"][:]
            coords_fold = f["coords/coords_center"][:]
            chroms_fold = f["coords/coords_chrom"][:]
            if (chroms is None) or (coords is None):
                chroms = chroms_fold
                coords = coords_fold
            else:
                assert(np.array_equal(chroms, chroms_fold))
                assert(np.array_equal(coords, coords_fold))
        data_ss.append(logcts[indicators == 1])
    
    mat_ss = np.column_stack(data_ss)

    data_xs = []
    for p in hdf5_paths_xs:    
        with h5py.File(p) as f:
            indicators = f["coords/coords_peak"][:]
            logcts = f["predictions/logcounts"][:]
            coords_fold = f["coords/coords_center"][:]
            chroms_fold = f["coords/coords_chrom"][:]
            if (chroms is None) or (coords is None):
                chroms = chroms_fold
                coords = coords_fold
            else:
                assert(np.array_equal(chroms, chroms_fold))
                assert(np.array_equal(coords, coords_fold))
        data_xs.append(logcts[indicators == 1])
    
    mat_xs = np.column_stack(data_xs)

    chroms = chroms.astype(str)

    true_counts_data = {}
    with open(true_counts_path) as f:
        next(f)
        for line in f:
            chrom, summit, logcts = line.rstrip("\n").split("\t")
            summit = int(summit)
            logcts = float(logcts)

            true_counts_data[(chrom, summit)] = logcts

    _, chrom_idx = np.unique(chroms, return_index=True)
    chrom_order = chroms[np.sort(chrom_idx)]
    counts_data = {}
    for c in chrom_order:
        inds_chrom, = np.nonzero(chroms == c)
        coords_chrom_unsorted = coords[inds_chrom]
        sort_order = np.argsort(coords_chrom_unsorted)
        coords_chrom = coords_chrom_unsorted[sort_order]
        # print(inds_chrom)####
        # print(coords_chrom) ####
        # print(sort_order) ####
        inds = inds_chrom[sort_order]

        # mat_chrom = mat[inds,:]
        logcounts_true = np.fromiter((true_counts_data[(c, i)] for i in coords_chrom), float, count=np.size(inds))
        record = {
            "center": coords_chrom,
            "logcounts_ss": mat_ss[inds,:],
            "logcounts_xs": mat_xs[inds,:],
            "logcounts_true": logcounts_true
        }
        counts_data[c] = record

    return counts_data

def load_gtf_tss(gtf_path, assembly):
    transcripts = {}
    with gzip.open(gtf_path, "rt") as i:
        for line in i:
            if line.startswith("#"):
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = line.rstrip('\n').split('\t')

            if feature != "transcript":
                continue

            if assembly == "mm10":
                chrom = f"chr{chrom}"

            if strand == "+":
                pos = int(start) - 1
            elif strand == "-":
                pos = int(end) - 1
            else:
                continue
            
            attributes = [i.strip().split(" ", 1) for i in attributes.rstrip(";").split(";")]
            # print([i for i in attributes if len(i) > 2])####
            attributes = {k: v.strip("\"") for k, v in attributes}            
            gene_id = attributes["gene_id"]
            gene_name = attributes.get("gene_name", gene_id)
            transcript_id = attributes["transcript_id"]
            transcript_name = attributes.get("transcript_name", transcript_id)

            record = {
                "gene_name": gene_name, 
                "gene_id": gene_id,
                "transcript_name": transcript_name,
                "transcript_id": transcript_id,
                "chrom": chrom,
                "pos": pos,
            }
            transcripts.setdefault(chrom, []).append(record)

    return transcripts

def get_primary_tss(transcripts, counts_data):
    primary_tss_data = {}
    for c, transcripts_chrom in transcripts.items():
        primary_tss_data.setdefault(c, {})
        counts_chrom = counts_data[c]
        if counts_chrom not in counts_data:
            continue
        for t in transcripts_chrom:
            pos_min = t["pos"] - 500
            pos_max = t["pos"] + 500

            idx_min = np.searchsorted(counts_chrom["center"], pos_min, side="left")
            idx_max = np.searchsorted(counts_chrom["center"], pos_max, side="right") 

            if idx_max <= idx_min:
                continue

            # tss_peak_idxs = counts_chrom["center"][idx_min:idx_max]
            best_tss_value = np.max(counts_chrom["logcounts_true"][idx_min:idx_max])
            record = t.copy()
            record["tss_count"] = best_tss_value
            
            gid = t["gene_id"]
            if gid not in primary_tss_data[c]:
                primary_tss_data[c][gid] = record
            elif primary_tss_data[c][gid]["tss_count"] < best_tss_value:
                primary_tss_data[c][gid] = record

    return primary_tss_data

def compute_scores(primary_tss_data, counts_data, slop, gamma):
    score_data = {}
    for chrom, t in primary_tss_data.items():
        score_data[chrom] = {}
        c = counts_data[chrom]
        genes = list(t.keys())
        # scores = np.zeros(len(genes), 5)
        for g in genes:
            tss_pos = t[g]["pos"]
            pos_min = tss_pos - slop
            pos_max = tss_pos + slop
            idx_min = np.searchsorted(c["center"], pos_min, side="left")
            idx_max = np.searchsorted(c["center"], pos_max, side="right") 

            if idx_max <= idx_min:
                ss_scores = np.full(5, 0.)
                xs_scores = np.full(5, 0.)
                score_data[chrom][g] = (ss_scores, xs_scores)
                continue

            ss_preds = c["logcounts_ss"][idx_min:idx_max,:]
            xs_preds = c["logcounts_xs"][idx_min:idx_max,:]
            peak_pos = c["pos"][idx_min:idx_max]

            dists = np.abs(peak_pos - tss_pos)
            weights = np.log(get_powerlaw_at_distance(dists, gamma))[:, np.newaxis]

            ss_scores = logsumexp(ss_preds + weights, axis=0)
            xs_scores = logsumexp(xs_preds + weights, axis=0)

            score_data[chrom][g] = (ss_scores, xs_scores)

    return score_data

def main(ss_paths, xs_paths, true_counts_path, gtf_path, out_path, assembly, slop, gamma):
    h5_ss_paths = [os.path.join(i, "out_predictions.h5") for i in ss_paths]
    h5_xs_paths = [os.path.join(i, "out_predictions.h5") for i in xs_paths]

    transcripts = load_gtf_tss(gtf_path, assembly)
    counts_data = load_logcts(h5_ss_paths, h5_xs_paths, true_counts_path)
    primary_tss_data = get_primary_tss(transcripts, counts_data)
    score_data = compute_scores(primary_tss_data, counts_data, slop, gamma)

    chroms = list(counts_data.keys())
    records = []
    scores_ss_all = []
    scores_xs_all = []
    for c in chroms:
        scores = score_data[c]
        for gene, tss_data in primary_tss_data[c].items():
            record = tss_data.copy()
            scores_ss, scores_xs = scores[gene]
            ss_data = {f"score_ss_{i}": v for i, v in enumerate(scores_ss)}
            xs_data = {f"score_ss_{i}": v for i, v in enumerate(scores_xs)}
            record.update(ss_data)
            record.update(xs_data)

            records.append(record)
            scores_ss_all.append(scores_ss)
            scores_xs_all.append(scores_xs)

    scores_df = pd.DataFrame.from_records(records)
    
    scores_mean_ss = np.mean(scores_ss_all, axis=1)
    scores_mean_xs = np.mean(scores_xs_all, axis=1)

    ss_n = scores_ss_all.shape[1]
    xs_n = scores_xs_all.shape[1]

    scores_ss_qn = rankdata(scores_ss_all, axis=0)
    scores_xs_qn = rankdata(scores_xs_all, axis=0)

    u_ss_qn, pval_qn = mannwhitneyu(scores_ss_qn, scores_xs_qn, method="exact", axis=1)
    u_xs_qn = ss_n * xs_n - u_ss_qn 
    nlp_qn = -np.log10(pval_qn)

    _, qval_qn = fdrcorrection(pval_qn)
    nlq_qn = -np.log10(qval_qn)

    additional_cols = {
        "score_ss_mean": scores_mean_ss,
        "score_xs_mean": scores_mean_xs,
        "u_stat_ss": u_ss_qn,
        "u_stat_xs": u_xs_qn,
        "-log10p": nlp_qn,
        "-log10q": nlq_qn
    }
    df_additional = pd.DataFrame(additional_cols)
    scores_df = pd.concat([scores_df, df_additional], axis=1)
    scores_df.sort_values("-log10q", inplace=True, ascending=False)

    scores_df.to_csv(out_path, sep="\t", index=False)

out_path, = snakemake.output

ss_paths = snakemake.input["ss_data"]
xs_paths = snakemake.input["xs_data"]
true_counts_path = snakemake.input["true_counts"]
gtf_path = snakemake.input["gtf"]

assembly = snakemake.wildcards["assembly"]

slop = int(snakemake.params["slop"])
gamma = snakemake.params["gamma"]

log_path = snakemake.log["score"]

main(ss_paths, xs_paths, true_counts_path, gtf_path, out_path, assembly, slop, gamma)