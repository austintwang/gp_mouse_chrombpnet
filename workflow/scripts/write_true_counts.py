from types import SimpleNamespace
import os

import numpy as np

import chrombpnet.training.data_generators.initializers as initializers


def output_counts(test_generator):
    num_batches=len(test_generator)
    true_counts_sum = []
    coordinates = []

    for idx in range(num_batches):
        if idx%100==0:
            print(str(idx)+'/'+str(num_batches))
        
        X,y,coords=test_generator[idx]

        true_counts_sum.extend(y[1][:,0])
        coordinates.extend(coords)

    return np.array(true_counts_sum), np.array(coordinates)

def main(model_dir, genome_path, bw_path, peaks_path, folds_path, out_path):
    model_path = os.path.join(model_dir, "chrombpnet_wo_bias.h5")
    args = SimpleNamespace(
        genome = genome_path,
        bigwig = bw_path,
        peaks = peaks_path,
        nonpeaks = "None",
        output_prefix = "tmp/",
        chr_fold_path = folds_path,
        trackables = ['loss','val_loss'],
        model_h5 = model_path,
        batch_size = 512,
        seed = 1234,
        inputlen = 2114,
        outputlen = 1000
    )
    test_generator = initializers.initialize_generators(args, mode="test", parameters=None, return_coords=True)
    true_counts_sum, coords = output_counts(test_generator)

    print(true_counts_sum) ####
    print(coords) ####

    num_examples=len(coords)

    coords_chrom_dset =  [str(coords[i][0]) for i in range(num_examples)]
    coords_center_dset =  [int(coords[i][1]) for i in range(num_examples)]
    coords_peak_dset =  [int(coords[i][3]) for i in range(num_examples)]
    true_counts_dset = [int(true_counts_sum[i]) for i in range(num_examples)]

    with open(out_path, "w") as f:
        f.write("chrom\tsummit_pos\tlog1p_true_counts\n")
        for a, b, c, d in zip(coords_chrom_dset, coords_center_dset, coords_peak_dset, true_counts_dset):
            if c == '1':
                f.write(f"{a}\t{b}\t{d}\n")


out_path, = snakemake.output

model_dir = snakemake.input["model_dir"]
genome_path = snakemake.input["genome"]
bw_path = snakemake.input["bigwig"]
peaks_path = snakemake.input["peaks"]
folds_path = snakemake.input["folds"]

main(model_dir, genome_path, bw_path, peaks_path, folds_path, out_path)