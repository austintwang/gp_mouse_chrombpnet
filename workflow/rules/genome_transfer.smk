def other_asm(assembly):
    if assembly == "mm10":
        return "cavpor3"
    elif assembly == "cavpor3":
        return "mm10"
    else:
        raise ValueError

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

rule filter_edge_peaks:
    """
    Build dummy test fold containing all peaks
    """
    input:
        peaks = "inputs/assembly/{assembly}/clusters/{cluster}/peaks.narrowPeak",
        bigwig = "inputs/assembly/{assembly}/clusters/{cluster}/coverage.bw",
    output:
        "results/assembly/{assembly}/clusters/{cluster}/peaks_edge_filtered.bed"
    conda:
        "../envs/chrombpnet.yaml"
    script:
        "../scripts/filter_edge_peaks.py"

rule peaks_global_chrombpnet_predict:
    """
    Run chrombpnet on peak regions
    """
    input:
        fasta = "inputs/assembly/{assembly}/genome.fa",
        peaks = "results/assembly/{assembly}/clusters/{cluster}/peaks_edge_filtered.bed",
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
        peaks = "results/assembly/{assembly}/clusters/{cluster}/peaks_edge_filtered.bed",
        train_dir = lambda w: f"results/assembly/{other_asm(w.assembly)}/clusters/{w.cluster}/folds/{w.fold}/train",
        bigwig = "inputs/assembly/{assembly}/clusters/{cluster}/coverage.bw",
        folds = "genomes/{assembly}_folds_all.json"
    output:
        directory("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/xs")
    log:
        transfer = "logs/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/xs.txt"

rule score_peaks:
    """
    Score peaks across species using genome transfer
    """
    input:
        ss_data = expand("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/ss", fold=config["folds_used"], allow_missing=True),
        xs_data = expand("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/xs", fold=config["folds_used"], allow_missing=True)
    output:
        "results/assembly/{assembly}/clusters/{cluster}/transfer/peak_scores.tsv"
    log:
        score = "logs/assembly/{assembly}/clusters/{cluster}/transfer/peak_scores.log"
    conda:
        "../envs/genome_transfer.yaml"
    script:
        "../scripts/score_peaks.py"

rule score_peaks_u:
    """
    Score peaks across species using genome transfer (Mann-Whitney U test)
    """
    input:
        ss_data = expand("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/ss", fold=config["folds_used"], allow_missing=True),
        xs_data = expand("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/transfer/xs", fold=config["folds_used"], allow_missing=True)
    output:
        "results/assembly/{assembly}/clusters/{cluster}/transfer/peak_scores_u.tsv"
    log:
        score = "logs/assembly/{assembly}/clusters/{cluster}/transfer/peak_scores_u.log"
    conda:
        "../envs/genome_transfer.yaml"
    script:
        "../scripts/score_peaks.py"

rule extract_true_counts:
    """
    Extract true count data for peaks
    """
    input:
        model_dir = "results/assembly/{assembly}/clusters/{cluster}/folds/0/train",
        genome = "inputs/assembly/{assembly}/genome.fa",
        peaks = "results/assembly/{assembly}/clusters/{cluster}/peaks_edge_filtered.bed",
        bigwig = "inputs/assembly/{assembly}/clusters/{cluster}/coverage.bw",
        folds = "genomes/{assembly}_folds_all.json"
    output:
        "results/assembly/{assembly}/clusters/{cluster}/transfer/true_counts.tsv"
    conda:
        "../envs/chrombpnet.yaml"
    resources:
        mem_mb = 40000
    script:
        "../scripts/write_true_counts.py"

rule merge_peak_scores:
    """
    Merge peak scoring data
    """
    input:
        scores = expand("results/assembly/{assembly}/clusters/{cluster}/transfer/peak_scores.tsv", cluster=clusters_l3, allow_missing=True),
        counts = expand("results/assembly/{assembly}/clusters/{cluster}/transfer/true_counts.tsv", cluster=clusters_l3, allow_missing=True)
    output:
        data = "results/assembly/{assembly}/scores_merged/data.tsv",
        summary = "results/assembly/{assembly}/scores_merged/summary.tsv"
    params:
        clusters = clusters_l3
    conda:
        "../envs/genome_transfer.yaml"
    resources:
        mem_mb = 40000
    script:
        "../scripts/merge_peak_scores.py"

rule merge_peak_scores_u:
    """
    Merge peak scoring data (U test)
    """
    input:
        scores = expand("results/assembly/{assembly}/clusters/{cluster}/transfer/peak_scores_u.tsv", cluster=clusters_l3, allow_missing=True),
        counts = expand("results/assembly/{assembly}/clusters/{cluster}/transfer/true_counts.tsv", cluster=clusters_l3, allow_missing=True)
    output:
        data = "results/assembly/{assembly}/scores_merged/data_u.tsv",
        summary = "results/assembly/{assembly}/scores_merged/summary_u.tsv"
    params:
        clusters = clusters_l3
    conda:
        "../envs/genome_transfer.yaml"
    resources:
        mem_mb = 40000
    script:
        "../scripts/merge_peak_scores_u.py"