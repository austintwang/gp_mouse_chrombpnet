rule build_folds_all:
    """
    Build dummy test fold containing all peaks
    """
    input:
        "inputs/assembly/{assembly}/folds"
    output:
        "genomes/{assembly}_folds_all.json"
    conda:
        "../envs/fetch.yaml"
    script:
        "../scripts/build_folds_all.py"

rule peaks_global_chrombpnet_predict:
    """
    Run chrombpnet on peak regions
    """
    input:
        fasta = "inputs/assembly/{assembly}/genome.fa",
        peaks = "inputs/assembly/{assembly}/clusters/{cluster}/peaks.narrowPeak",
        train_dir = "results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/train",
        bigwig = "inputs/assembly/{assembly}/clusters/{cluster}/coverage.bw",
        folds = "genomes/{assembly}_folds_all.json"
    output:
        directory("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/ss")
    log:
        transfer = "logs/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/ss.txt"
    conda:
        "../envs/chrombpnet.yaml"
    resources:
        mem_mb = 40000,
        runtime = "24:00:00",
        partitions = "akundaje,gpu,owners",
        gpu_cmd = "-G 1 -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'"
    shell:
        "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}; "
        "mkdir -p {output}; "
        "chrombpnet_predict "
        "--genome={input.fasta} "
        "--bigwig={input.bigwig} "
        "--peaks={input.peaks} " 
        "--chr_fold_path={input.folds} "
        "--inputlen=2114 "
        "--outputlen=1000 "
        "--output_prefix={output}/out "
        "--batch_size=256 "
        "--model_h5={input.train_dir}/chrombpnet_wo_bias.h5 &> {log.transfer}"
        
use rule peaks_global_chrombpnet_predict as peaks_global_chrombpnet_predict_xs with:
    input:
        fasta = "inputs/assembly/{assembly}/genome.fa",
        peaks = "inputs/assembly/{assembly}/clusters/{cluster}/peaks.narrowPeak",
        model_dir = lambda w: f"results/assembly/{other_asm(w.assembly)}/clusters/{w.cluster}/folds/{w.fold}/train",
        bigwig = "inputs/assembly/{assembly}/clusters/{cluster}/coverage.bw",
        folds = "genomes/{assembly}_folds_all.json"
    output:
        directory("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/xs")
    log:
        transfer = "logs/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/xs.txt"

rule xs_global_visualize:
    """
    Visualize projection performance
    """
    input:
        ss_data = "results_merged/{assembly}/xs_chrombpnet_global/markers_{cluster}/fold_{fold}_peaks_chrombpnet_predict_ss",
        xs_data = "results_merged/{assembly}/xs_chrombpnet_global/markers_{cluster}/fold_{fold}_peaks_chrombpnet_predict_xs"
    output:
        directory("results_merged/{assembly}/xs_chrombpnet_global/markers_{cluster}/fold_{fold}_peaks_chrombpnet_compare")
    conda:
        "../envs/chrombpnet_legacy.yaml"
    script:
        "../scripts/xs_chrombpnet_global_visualize.py"