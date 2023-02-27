import os
import itertools
import numpy as np
import h5py
# from scipy.stats import t
from scipy.stats import mannwhitneyu, rankdata
from scipy.spatial import distance
from scipy.special import comb
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd

def load_profiles(hdf5_paths_ss, hdf5_paths_xs):
    chroms = None
    coords = None
    profs_ss = []
    profs_xs = []

    for p in hdf5_paths_ss:
        with h5py.File(p) as f:
            indicators = f["coords/coords_peak"][:]
            coords_fold = f["coords/coords_center"][:]
            chroms_fold = f["coords/coords_chrom"][:]
            profs = f["predictions/profs"][:]
            if (chroms is None) or (coords is None):
                chroms = chroms_fold
                coords = coords_fold
            else:
                assert(np.array_equal(chroms, chroms_fold))
                assert(np.array_equal(coords, coords_fold))
            profs = profs[indicators==1,:]
            assert(len(profs) == len(coords))
        profs_ss.append(profs)

    for p in hdf5_paths_xs:
        with h5py.File(p) as f:
            indicators = f["coords/coords_peak"][:]
            coords_fold = f["coords/coords_center"][:]
            chroms_fold = f["coords/coords_chrom"][:]
            profs = f["predictions/profs"][:]
            if (chroms is None) or (coords is None):
                chroms = chroms_fold
                coords = coords_fold
            else:
                assert(np.array_equal(chroms, chroms_fold))
                assert(np.array_equal(coords, coords_fold))
            profs = profs[indicators==1,:]
            assert(len(profs) == len(coords))
        profs_xs.append(profs)

    return profs_ss, profs_xs, coords, chroms

def pairwise_jsds(profs_ss, profs_xs):
    num_peaks = len(profs_ss[0])

    # n = len(profs_ss)
    # ss_self_jsds = np.zeros((n*(n-1)/2, num_peaks),)
    # for i, pair in enumerate(itertools.combinations(profs_ss, 2)):
    #     a, b = pair
    #     jsd = distance.jensenshannon(a, b, axis=1)
    #     ss_self_jsds[i,:] = jsd

    # n = len(profs_xs)
    # xs_self_jsds = np.zeros((n*(n-1)/2, num_peaks),)
    # for i, pair in enumerate(itertools.combinations(profs_xs, 2)):
    #     a, b = pair
    #     jsd = distance.jensenshannon(a, b, axis=1)
    #     xs_self_jsds[i,:] = jsd

    n = len(profs_ss)
    ss_self_jsds = np.zeros((num_peaks,n,n),)
    for i, a in enumerate(profs_ss):
        for j, b in enumerate(profs_ss):
            jsd = distance.jensenshannon(a, b, axis=1)
            ss_self_jsds[:,i,j] = jsd

    n = len(profs_xs)
    xs_self_jsds = np.zeros((num_peaks,n,n),)
    for i, a in enumerate(profs_xs):
        for j, b in enumerate(profs_xs):
            jsd = distance.jensenshannon(a, b, axis=1)
            xs_self_jsds[:,i,j] = jsd

    cross_jsds = np.zeros((num_peaks,len(profs_ss),len(profs_xs)),)
    for i, a in enumerate(profs_ss):
        for j, b in enumerate(profs_xs):
            jsd = distance.jensenshannon(a, b, axis=1)
            cross_jsds[:,i,j] = jsd

    return ss_self_jsds, xs_self_jsds, cross_jsds

def permutation_test_energy(ss_self_jsds, xs_self_jsds, cross_jsds, energy_scores):
    num_ss, num_xs = cross_jsds.shape[1:3]
    num_total = num_ss + num_xs
    jsds_all = np.block([
        [ss_self_jsds, cross_jsds],
        [cross_jsds.transpose(0,2,1), xs_self_jsds],
    ])
    num_shuffles = comb(num_total, num_ss, exact=True)
    shuffles = np.full((1,num_total,num_shuffles), -1)
    for i, c in enumerate(itertools.combinations(range(num_total), num_ss)):
        shuffles[:,c,i] = 1
    
    scores = (shuffles * (jsds_all @ shuffles)).sum(axis=1)
    scores_sorted = np.sort(scores, axis=1)
    nulls = scores_sorted - energy_scores
    pos = np.apply_along_axis(np.searchsorted, 1, nulls, 0)
    pvals = 1 - pos / num_shuffles

    return pvals

def main(ss_paths, xs_paths, out_path):
    h5_ss_paths = [os.path.join(i, "out_predictions.h5") for i in ss_paths]
    h5_xs_paths = [os.path.join(i, "out_predictions.h5") for i in xs_paths]

    profs_ss, profs_xs, coords, chroms = load_profiles(h5_ss_paths, h5_xs_paths)

    ss_self_jsds, xs_self_jsds, cross_jsds = pairwise_jsds(profs_ss, profs_xs)
    ss_self_jsd_mean = np.mean(ss_self_jsds, axis=(1,2))
    xs_self_jsd_mean = np.mean(xs_self_jsds, axis=(1,2))
    cross_jsd_mean = np.mean(cross_jsds, axis=(1,2))
    prof_e_dist = 2 * cross_jsd_mean - ss_self_jsd_mean - xs_self_jsd_mean

    pval = permutation_test_energy(ss_self_jsds, xs_self_jsds, cross_jsds, prof_e_dist)
    nlp = -np.log10(pval)

    _, qval = fdrcorrection(pval)
    nlq = -np.log10(qval)
   
    cols = {
        "chrom": chroms, 
        "summit_pos": coords, 
        "ss_self_jsd_mean": ss_self_jsd_mean,
        "xs_self_jsd_mean": xs_self_jsd_mean,
        "cross_jsd_mean": cross_jsd_mean,
        "prof_e_dist": prof_e_dist,
        "prof_nlp": nlp,
        "prof_nlq": nlq,
    }
    # for k, v in cols.items(): ####
    #     print(k)
    #     print(v.shape)
    #     print(v) 
    table = pd.DataFrame(cols)
    table.sort_values("prof_e_dist", inplace=True, ascending=False)

    table.to_csv(out_path, sep="\t", index=False)


out_path, = snakemake.output

ss_paths = snakemake.input["ss_data"]
xs_paths = snakemake.input["xs_data"]

main(ss_paths, xs_paths, out_path)