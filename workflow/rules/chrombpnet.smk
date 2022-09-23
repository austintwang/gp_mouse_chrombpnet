
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
        mem_mb = 30000,
        runtime = "24:00",
        gpu_cmd = "-G 1 -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'"
    conda:
        "../envs/chrombpnet.yaml"
    shell:
        "export CUDNN=cudnn-8.1_cuda11.2; "
        "export cuda=cuda-11.2; "
        "export LD_LIBRARY_PATH=/usr/local/$cuda/lib64:/usr/local/$CUDNN/lib64:/usr/local/$CUDNN/include:/usr/local/$cuda/extras/CUPTI/lib64:/usr/local/lib:$LD_LIBRARY_PATH; "
        "export PATH=/usr/local/$cuda/bin:$PATH; "
        "export CUDA_HOME=/usr/local/$cuda; "
        "export CPATH='/usr/local/$CUDNN/include:${{CPATH}}'; "
        "export LIBRARY_PATH='/usr/local/$CUDNN/lib64:${{LIBRARY_PATH}}'; "
        "export CPLUS_INCLUDE_PATH=/usr/local/$cuda/include; "
        "mkdir -p {output}; "
        "step6_train_chrombpnet_model.sh {input.fasta} {input.bigwig} {input.peaks} {input.negatives} " 
        "{input.folds}/fold_{wildcards.fold}.json {input.bias_model} {output} ATAC_PE {log.train}"
