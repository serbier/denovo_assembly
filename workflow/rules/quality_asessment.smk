rule busco:
	input:
		scaffolded_draft = 'results/{sample}/scaffolding/{sample}_scaffolded.fasta'
	output:
		directory('results/{sample}/BUSCO/')
	threads:
		10
	conda:
		"denovo"
	shell:
		"""
		{config[BUSCO]} -i {input.scaffolded_draft} \
				-m geno \
				-o {output} \
				-c {threads} \
				--offline
		"""