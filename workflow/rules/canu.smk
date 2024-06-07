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
        "denovo_assembly_pipeline"
    shell: 
        """
        canu -trim-assemble \
        -p {wildcards.sample} \
        -d {output} \
        genomeSize={params.genome_size}m \
        useGrid=false \
        maxThreads={threads} \
        -corrected -nanopore {input} 1> {log}
        """

rule NGSEP_denovo_assembly:
    input:
        corrected_reads = get_corrected_reads
    output:
        'results/{sample}/assembly/NGSEP/{sample}.fa'
    log:
        'results/{sample}/assembly/NGSEP/{sample}_Assembler.log'
    threads:
        50
    conda:
        "denovo"
    shell:
        """
        java -Xmx300G -jar {config[NGSEP]} Assembler \
            -i {input} \
            -o results/{wildcards.sample}/assembly/NGSEP/{wildcards.sample} \
            -m 3000 \
            -t {threads} 2> {log}
        """