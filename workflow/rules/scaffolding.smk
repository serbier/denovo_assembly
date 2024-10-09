rule index_curated_draft:
	input:
		curated_draft = 'results/{sample}/purge/{sample}.fasta'
	output:
		indexes = multiext('results/{sample}/purge/{sample}.fasta',
				".amb",
				".ann",
				".bwt",
				".pac",
				".sa")
	shell:
		"""
		bwa index {input.curated_draft}
		"""

rule map_markers_to_assembly:
	input:
		markers_sequence = config['markers_scaffolding'],
		curated_draft = 'results/{sample}/purge/{sample}.fasta',
		indexes = multiext('results/{sample}/purge/{sample}.fasta',
				".amb",
				".ann",
				".bwt",
				".pac",
				".sa")
	output:
		mapping = 'results/{sample}/scaffolding/{sample}_SR_F2.sam'
	conda:
		"denovo_assembly_pipeline"
	threads:
		10
	log:
		'results/{sample}/scaffolding/{sample}_SR_F2.log'
	shell:
		"""
		bwa aln -t {threads} {input.curated_draft} {input.markers_sequence} > \
			results/{wildcards.sample}/scaffolding/{wildcards.sample}_SR_F2.sai 2> \
			{log} && \
		bwa samse {input.curated_draft} results/{wildcards.sample}/scaffolding/{wildcards.sample}_SR_F2.sai \
			{input.markers_sequence} \
			> {output.mapping} && rm results/{wildcards.sample}/scaffolding/{wildcards.sample}_SR_F2.sai
		"""
rule get_AGP_from_curated_draft:
	input:
		curated_draft = 'results/{sample}/purge/{sample}.fasta',
	output:
		agp = 'results/{sample}/scaffolding/{sample}.agp'
	shell:
		"""
		{config[abyss-fatoagp]} {input.curated_draft} | \
			sed 's/^scaffold//g'| \
			sed 's/contigC/C/g' | \
			sed 's/contigS/S/g' | \
			sed 's/_0//g'> {output.agp}
		"""

rule run_chromonomer:
	input:
		linkage_map = config['linkage_map'],
		mapping = 'results/{sample}/scaffolding/{sample}_SR_F2.sam',
		agp = 'results/{sample}/scaffolding/{sample}.agp'
	output:
		outdir = directory('results/{sample}/scaffolding/chromonomer'),
		agp = 'results/{sample}/scaffolding/chromonomer/CHRR_genome.agp'
	log:
		'results/{sample}/scaffolding/chromonomer/{sample}_chromonomer.log'
	shell:
		"""
		{config[chromonomer]} -p {input.linkage_map} \
			-s {input.mapping} \
			-a {input.agp} \
			--data_version SR_F2 \
			-o {output.outdir} 2> {log}
		"""

rule get_scaffolded_assembly:
	input:
		agp = 'results/{sample}/scaffolding/chromonomer/CHRR_genome.agp',
		curated_draft = 'results/{sample}/purge/{sample}.fasta',
	output:
		'results/{sample}/scaffolding/{sample}_scaffolded.fasta'
	conda:
		"ragtag"
	log:
		'results/{sample}/scaffolding/{sample}_agp2fasta.log'

	shell:
		"""
		ragtag_agp2fa.py {input.agp} {input.curated_draft} > {output} 2> {log}
		"""

rule ref_based_scaffolding:
	input:
		ref = config['scaffolding']['ref'],
		curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta",
	output:
		multiext(f"{base_dir}/scaffolding/ragtag/{{sample}}/ragtag.scaffold",".agp", ".fasta", '.stats')
		
	conda:
		"../envs/ragtag.yaml"
	log:
		'logs/{sample}/scaffolding/ragtag_scaff.log'
	threads:
		30
	shell:
		"""
		ragtag.py scaffold {input.ref} {input.curated_draft} -o {base_dir}/scaffolding/ragtag/{wildcards.sample} -t {threads} 2> {log}
		"""



