""" If the sample was sequenced with multiple runs 
is necessary to merge the reads into a single file 
otherwise use the basecalled file """

rule merge_raw_reads:
    input:
        get_read_list_per_sample
    output:
        f"{base_dir}/basecalling/{{sample}}/merged.fastq.gz"
    conda:
        "ngs"
    shell:
        """
        zcat {input}| bgzip > {output}
        """
rule bring_raw_reads:
    input:
        reads = get_raw_reads
    output:
        reads = temp("data/raw_reads/{sample}.fastq.gz"),
        index = temp("data/raw_reads/{sample}.fastq.gz.fai"),
    conda:
        "ngs"
    shell:
        """
        cp {input.reads} {output.reads} && \
        samtools faidx {output.reads}
        """