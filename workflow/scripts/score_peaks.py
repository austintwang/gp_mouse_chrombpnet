import os
import numpy as np
import h5py
from scipy.stats import t
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd

def load_logcts(hdf5_paths):
    chroms = None
    coords = None
    data = []
    for p in hdf5_paths:    
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
        data.append(logcts[indicators == 1])
    
    mat = np.column_stack(data)

    return mat, chroms.astype(str), coords

def main(ss_paths, xs_paths, out_path, log_path):
    h5_ss_paths = [os.path.join(i, "out_predictions.h5") for i in ss_paths]
    h5_xs_paths = [os.path.join(i, "out_predictions.h5") for i in xs_paths]

    ss_data, ss_chroms, ss_coords = load_logcts(h5_ss_paths)
    xs_data, xs_chroms, xs_coords = load_logcts(h5_xs_paths)

    assert(np.array_equal(ss_chroms, xs_chroms))
    assert(np.array_equal(ss_coords, xs_coords))

    chroms = ss_chroms
    coords = ss_coords

    ss_data -= ss_data.mean()
    xs_data -= xs_data.mean()

    ss_n = ss_data.shape[1]
    xs_n = xs_data.shape[1]

    ss_mean = np.mean(ss_data, axis=1)
    xs_mean = np.mean(xs_data, axis=1)

    ss_std = np.std(ss_data, axis=1, ddof=1)
    xs_std = np.std(xs_data, axis=1, ddof=1)

    t_denom = np.sqrt(ss_std**2/ss_n + xs_std**2/xs_n)
    t_stat = (xs_mean - ss_mean) / t_denom

    df_num = (ss_std**2/ss_n + xs_std**2/xs_n)**2
    df_denom = ((ss_std**2/ss_n)**2/(ss_n-1) + (xs_std**2/xs_n)**2/(xs_n-1))
    df = df_num / df_denom

    nlp = -(t.logsf(np.abs(t_stat), df) + np.log(2)) / np.log(10)
    pval = 10**(-nlp)
    _, qval = fdrcorrection(pval)
    nlq = -np.log10(qval)

    cols = {
        "chrom": chroms, 
        "summit_pos": coords, 
        "ss_mean": ss_mean,
        "xs_mean": xs_mean,
        "ss_std": ss_std,
        "xs_std": xs_std,
        "t_stat": t_stat,
        "est_df": df,
        "-log10p": nlp,
        "-log10q": nlq
    }
    # for k, v in cols.items(): ####
    #     print(k)
    #     print(v.shape)
    #     print(v) 
    table = pd.DataFrame(cols)
    table.sort_values("-log10q", inplace=True)

    table.to_csv(out_path, sep="\t", index=False)

    f10 = nlq >= 1
    f1 = nlq >= 2
    with open(log_path, "w") as f:
        f.write(f"10% FDR:\t{f10.sum()}\t{f10.mean()}\n")
        f.write(f"1% FDR:\t{f1.sum()}\t{f1.mean()}\n")

out_path, = snakemake.output

ss_paths = snakemake.input["ss_data"]
xs_paths = snakemake.input["xs_data"]

log_path = snakemake.log["score"]

main(ss_paths, xs_paths, out_path, log_path)