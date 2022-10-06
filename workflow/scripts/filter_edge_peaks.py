import pyBigWig
import pandas as pd
from chrombpnet.helpers.hyperparameters import param_utils as param_utils

def main(peaks_path, bw_path, out_path):
    peaks =  pd.read_csv(peaks_path,
                           sep='\t',
                           header=None,
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
    bw = pyBigWig.open(bw_path) 

    peaks = param_utils.filter_edge_regions(peaks, bw, 2114+2*500, peaks_bool=1)

    peaks.to_csv(out_path, sep="\t",  header=False, index=False)

out_path, = snakemake.output

peaks_path = snakemake.input["peaks"]
bw_path = snakemake.input["bigwig"]

main(peaks_path, bw_path, out_path)