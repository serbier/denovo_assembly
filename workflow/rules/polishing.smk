rule get_overlaps:
	input:
		draft = 'results/{sample}/assembly/NGSEP/{sample}.fa',
		reads = get_corrected_reads 
	output:
		overlaps = temp('results/{sample}/polish/overlaps/{sample}.overlaps.sam')
	conda:
		'denovo_assembly_pipeline'
	threads:
		50
	shell:
		"""
		minimap2 -a -t {threads} {input.draft} {input.reads}  > {output.overlaps}
		"""



rule polish_genome:
	input:
		overlaps = 'results/{sample}/polish/overlaps/{sample}.overlaps.sam',
		reads = get_corrected_reads ,
		draft = 'results/{sample}/assembly/NGSEP/{sample}.fa',
	output:
		polish_draft = 'results/{sample}/polish/racon/{sample}.polished.fa'
	conda:
		"denovo"
	threads:
		50
	log: 
		'results/{sample}/polish/racon/{sample}.polishing.log'

	shell:
		"""
		{config[racon]} -t {threads} {input.reads} {input.overlaps} {input.draft} > {output.polish_draft} 2> {log}
		"""

rule index_polish_genome:
	input:
		polish_draft = 'results/{sample}/polish/racon/{sample}.polished.fa'
	output:
		idx = 'results/{sample}/polish/racon/{sample}.polished.fa.fai'
	conda:
		'denovo_assembly_pipeline'
	shell:
		"""
		samtools faidx {input.polish_draft}
		"""

rule samtools_faidx:
    input:
        "data/{sample}/{sample}_merged_reads.fastq.gz",
    output:
        temp("data/{sample}/{sample}_merged_reads.fastq.gz.fai")
    params:
        extra="",
    wrapper:
        "v3.12.1/bio/samtools/faidx"


rule bring_assembly:
	input:
		directory(f'{base_dir}/assembly/flye/{{sample}}')
	output:
		temp(directory('polishing/assembly_base/flye/{sample}'))
	shell:
		"""
		cp -r {input} {output}
		"""

rule polish_medaka:
	input:
		basecalls = "data/{sample}/{sample}_merged_reads.fastq.gz",
		draft = directory('polishing/assembly_base/flye/{sample}')
	output:
		temp(directory("polishing/medaka/{sample}"))
	threads:
		50
	log:
		'logs/{sample}/polishing/medaka.log'
	conda:
		"tf2"
	shell:
		"""
		medaka_consensus \
			-i {input.basecalls} \
			-d polishing/assembly_base/flye/{wildcards.sample}/assembly.fasta\
			-o {output} \
			-t {threads} \
			-m r1041_e82_400bps_hac_v4.3.0 2> {log}
		"""


