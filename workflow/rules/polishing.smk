rule polish_medaka:
	input:
		basecalls = "data/raw_reads/{sample}.fastq.gz",
		draft = "data/assembly/{sample}.fasta",
		index = "data/assembly/{sample}.fasta.fai"
	output:
		consensus = f"{base_dir}/polishing/medaka/{{sample}}/consensus.fasta",
		mappings = f"{base_dir}/polishing/medaka/{{sample}}/calls_to_draft.bam"
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
			-d {input.draft}\
			-o tmp_polishing/{wildcards.sample} \
			-t {threads} \
			-m r1041_e82_400bps_hac_v4.3.0 2> {log} && \
		cp -r tmp_polishing/{wildcards.sample} {base_dir}/polishing/medaka/ && \
		rm -r tmp_polishing/{wildcards.sample}
		"""


