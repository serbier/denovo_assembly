import os, sys

rule canu_denovo_assembly:
    input:
        corrected_reads = get_corrected_reads
    
    output:
        directory('results/{sample}/assembly/canu')
    
    threads: 
        resources['CANU']['threads']

    params:
        genome_size = config['organism']['genome_size_mb']

    log:
        "results/{sample}/assembly/canu.log"

    conda:
        "../envs/canu.yaml"
    shell: 
        """
        canu -trim-assemble \
        -p {wildcards.sample} \
        -d {output} \
        genomeSize={params.genome_size}m \
        useGrid=false \
        maxThreads={threads} \
        -corrected -nanopore {input} 2> {log}
        """

rule NGSEP_denovo_assembly:
    input:
        corrected_reads = get_corrected_reads
    output:
        'results/{sample}/assembly/NGSEP/{sample}.fa'
    log:
        'results/{sample}/assembly/NGSEP/{sample}_Assembler.log'
    threads:
        resources['NGSEP']['threads']
    conda:
        "../envs/NGSEP.yaml"
    shell:
        """
        java -Xmx300G -jar {config[NGSEP]} Assembler \
            -i {input} \
            -o results/{wildcards.sample}/assembly/NGSEP/{wildcards.sample} \
            -m 3000 \
            -t {threads} 2> {log}
        """

rule get_raw_reads:
    input:
        raw_reads = get_read_list_per_sample
    output:
        temp("data/{sample}/{sample}_merged_reads.fastq.gz")
    conda:
        'ngs'
    shell:
        """
        zcat {input.raw_reads} | bgzip > {output}
        """

rule flye_denovo_assembly:
    input:
        raw_reads = "data/{sample}/{sample}_merged_reads.fastq.gz"
    output:
        temp(directory('assembly/flye/{sample}')),

    log:
        'logs/{sample}/assembly/flye.log'
    conda:
        "../envs/flye.yaml"
    threads:
        resources['flye']['threads']
    shell:
        """
        flye --nano-raw {input.raw_reads} \
        --out-dir assembly/flye/{wildcards.sample} \
        --threads {threads} 2> {log}
        """
        
#rule move_assembly:
#    input:
#        directory('assembly/flye/{sample}')
#    output:
#        directory(f'{base_dir}/assembly/flye/{{sample}}')
#    shell:
#        """
#        cp -r {input} {output} && \
#        rm -r {input}
#        """
    