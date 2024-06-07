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

		