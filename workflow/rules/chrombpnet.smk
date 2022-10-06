
rule chrombpnet_train:
    """
    Train ChromBPNet model on cluster
    """
    input:
        fasta = "inputs/assembly/{assembly}/genome.fa",
        bigwig = "inputs/assembly/{assembly}/clusters/{cluster}/coverage.bw",
        peaks = "inputs/assembly/{assembly}/clusters/{cluster}/peaks.narrowPeak",
        negatives = "inputs/assembly/{assembly}/clusters/{cluster}/folds/{fold}/negatives.bed",
        folds = "inputs/assembly/{assembly}/folds",
        bias_model = "inputs/bias_model.h5"
    output:
        directory("results/assembly/{assembly}/clusters/{cluster}/folds/{fold}/train")
    log:
        train = "logs/assembly/{assembly}/clusters/{cluster}/folds/{fold}/train.txt"
    resources:
        mem_mb = 40000,
        runtime = "24:00:00",
        partitions = "akundaje,gpu,owners",
        gpu_cmd = "-G 1 -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'"
    conda:
        "../envs/chrombpnet.yaml"
    shell:
        "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}; "
        "mkdir -p {output}; "
        "step6_train_chrombpnet_model.sh {input.fasta} {input.bigwig} {input.peaks} {input.negatives} " 
        "{input.folds}/fold_{wildcards.fold}.json {input.bias_model} {output} ATAC_PE &> {log.train}"
