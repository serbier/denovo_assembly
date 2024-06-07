rule overlaps_for_purge: 
	input:
		polish_draft = 'results/{sample}/polish/racon/{sample}.polished.fa',
		reads = get_corrected_reads ,
	output:
		overlaps = temp('results/{sample}/purge/overlaps/{sample}.overlaps.sam')
	conda:
		'denovo_assembly_pipeline'
	threads:
		50

	shell:
		"""
		minimap2 -t {threads} -a -x map-ont \
		 {input.polish_draft} {input.reads} > {output.overlaps}
		"""

rule convert_overlaps_for_purge: 
	input:
		overlaps = 'results/{sample}/purge/overlaps/{sample}.overlaps.sam'
	output:
		overlaps = temp('results/{sample}/purge/overlaps/{sample}_sorted.overlaps.bam'),
		overlaps_idx = 'results/{sample}/purge/overlaps/{sample}_sorted.overlaps.bam.bai'
	conda:
		'denovo_assembly_pipeline'
	shell:
		"""
		samtools view -S -b {input.overlaps} > results/{wildcards.sample}/purge/overlaps/sample.bam
		samtools sort -o {output.overlaps} results/{wildcards.sample}/purge/overlaps/sample.bam
		rm results/{wildcards.sample}/purge/overlaps/sample.bam
		samtools index {output.overlaps}
		"""


rule get_read_histogram:
	input:
		overlaps = 'results/{sample}/purge/overlaps/{sample}_sorted.overlaps.bam',
		polish_draft = 'results/{sample}/polish/racon/{sample}.polished.fa',
		polish_draft_idx = 'results/{sample}/polish/racon/{sample}.polished.fa.fai'

	output:
		gen_coverage_stats =  'results/{sample}/purge/{sample}_polished.gencov',
		histogram = 'results/{sample}/purge/{sample}_polished_histogram.png',

	threads:
		50
	conda:
		'denovo_assembly_pipeline'

	shell:
		"""
		purge_haplotigs readhist -t {threads} \
			-b {input.overlaps} \
			-g {input.polish_draft} && \
			mv {wildcards.sample}_sorted.overlaps.bam.gencov {output.gen_coverage_stats} && \
			mv {wildcards.sample}_sorted.overlaps.bam.histogram.png {output.histogram}

		"""

rule coverage_stats_for_purge:
	input:
		gen_coverage_stats =  'results/{sample}/purge/{sample}_polished.gencov',
	output:
		contig_coverage =  'results/{sample}/purge/{sample}_polished_contig_coverage.csv',
	params:
		cutoffs = lambda wildcards: get_read_depth_cutoffs(wildcards),
	conda:
		'denovo_assembly_pipeline'
	log:
		'results/{sample}/purge/{sample}_polished_contig_coverage_stats.log'
	shell:
		"""
		purge_haplotigs contigcov -i {input.gen_coverage_stats} \
			-l {params[cutoffs][l]} \
			-m {params[cutoffs][m]} \
			-h {params[cutoffs][h]} && \
		mv coverage_stats.csv {output.contig_coverage}
		"""


rule purge_assembly:
	input:
		polish_draft = 'results/{sample}/polish/racon/{sample}.polished.fa',
		overlaps = 'results/{sample}/purge/overlaps/{sample}_sorted.overlaps.bam',
		overlaps_idx = 'results/{sample}/purge/overlaps/{sample}_sorted.overlaps.bam.bai',
		contig_coverage =  'results/{sample}/purge/{sample}_polished_contig_coverage.csv',
	output:
		dotplots = directory('results/{sample}/purge/dotplots'),
		curated_draft = 'results/{sample}/purge/{sample}.fasta',
		artefacts = 'results/{sample}/purge/{sample}.artefacts.fasta',
		log = 'results/{sample}/purge/{sample}.contig_associations.log',
		haplotigs = 'results/{sample}/purge/{sample}.haplotigs.fasta',
		reassignments = 'results/{sample}/purge/{sample}.reassignments.tsv'
	threads:
		30
	log:
		'results/{sample}/purge/{sample}_purge_haplotigs.log',
	conda:
		'denovo_assembly_pipeline'
	shadow:
		"full"
	shell:
		"""
		purge_haplotigs purge -g {input.polish_draft} \
			-t {threads} \
			-c {input.contig_coverage} \
			-o {wildcards.sample} \
			-d -b {input.overlaps} 2> {log} && \
		mkdir {output.dotplots} && \
		mv dotplots*/ {output.dotplots} && \
		mv tmp_purge_haplotigs results/{wildcards.sample}/purge && \
		mv {wildcards.sample}.fasta {output.curated_draft} && \
		mv {wildcards.sample}.artefacts.fasta {output.artefacts} && \
		mv {wildcards.sample}.contig_associations.log {output.log} && \
		mv {wildcards.sample}.haplotigs.fasta {output.haplotigs} && \
		mv {wildcards.sample}.reassignments.tsv {output.reassignments}
		"""