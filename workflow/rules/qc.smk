
#### QC ####
rule BamStats:
    """
    QC alignment statistics
    """
    input:
        bam = "results/alignments/{sample}.bam",
        idx = "results/alignments/{sample}.bam.bai"
    output:
        stats = "results/alignments/bamStats/{sample}.flagstat"
    log:
        "logs/BamStats/{sample}.log"
    wrapper:
        "0.70.0/bio/samtools/flagstat"

rule Coverage:
    """
    Calculate coverage with mosdepth
    """
    input:
        bam = "results/alignments/{sample}.bam",
        idx = "results/alignments/{sample}.bam.bai"
    output:
        "results/alignments/coverage/{sample}.mosdepth.summary.txt"
    log:
        "logs/Coverage/{sample}.log"
    conda:
        "../envs/depth.yaml"
    params:
        prefix = lambda w, output: output[0].split(os.extsep)[0],
        windowsize = 300
    threads:4
    shell: "mosdepth --threads {threads} --fast-mode --by {params.windowsize} --no-per-base {params.prefix} {input.bam}"