rule TE_annotation_scaffolded:
	input:
		purged_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta"
	output:
		te_annotation = f'{base_dir}/annotation/{{sample}}/TE/{{sample}}_TE_annotation_NGSEP.out'
	log:
		'logs/{sample}/annotation/NGSEP_TE.log'
	conda:
		"ngs"
	threads:
		30
	shell:
		"""
		java -Xmx50G  -jar {config[NGSEP]} TransposonsFinder \
			-i {input.purged_draft} \
			-o {output.te_annotation} \
			-d {config[TE_DB]} \
			-r 2 \
			-t {threads} 2> {log}
		"""

rule convert_NGSEP_TE_to_GFF:
	input:
		te_annotation = f'{base_dir}/annotation/{{sample}}/TE/{{sample}}_TE_annotation_NGSEP.out'
	output:
		te_annotation = f'{base_dir}/annotation/{{sample}}/TE/{{sample}}_TE_annotation_NGSEP.gff'
	shell:
		"""
		{config[NGSEP_TE2GFF_script]} {input.te_annotation} > \
			{output.te_annotation}
		"""

rule set_maker_config_files:
	input:
		draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta",
		te_annotation = f'{base_dir}/annotation/{{sample}}/TE/{{sample}}_TE_annotation_NGSEP.gff'


	output:
		multiext(f'{base_dir}/annotation/{{sample}}/maker/maker_',
		 'bopts.ctl',
		  'evm.ctl',
		    'exe.ctl',
			 'opts.ctl')
	shadow:
		'shallow'
	conda:
		"maker"
	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin
		maker -CTL && \
		awk -v est_path="{config[transcripts]}" \
			-v protein_path="{config[protein]}" \
			-v genome_path="{input.draft}" \
			-v te_path="{input.te_annotation}" \
			-v species="phaseouls_vulgaris_UI-111" -f {config[opts_maker_script]} \
			./maker_opts.ctl > ./tmp.ctl && \
		rm ./maker_opts.ctl && \
		mv ./maker_* {base_dir}/annotation/{wildcards.sample}/maker/ && \
		mv ./tmp.ctl {base_dir}/annotation/{wildcards.sample}/maker/maker_opts.ctl
		"""


rule run_maker:
	input:
		multiext(f'{base_dir}/annotation/{{sample}}/maker/maker_',
		 'bopts.ctl',
		  'evm.ctl',
		    'exe.ctl',
			 'opts.ctl')
	output:
		index = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.maker.output/{{sample}}_master_datastore_index.log'

	conda:
		"maker"
	threads:
		40
	log:
		'logs/{sample}/annotation/maker.log'
	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin && \
		{config[mpiexec]} -n {threads} maker -c 1 -q {base_dir}/annotation/{wildcards.sample}/maker/maker_opts.ctl \
		 {base_dir}/annotation/{wildcards.sample}/maker/maker_bopts.ctl \
		 {base_dir}/annotation/{wildcards.sample}/maker/maker_exe.ctl 2> {log} && \
		 mv ./{wildcards.sample}.maker.output {base_dir}/annotation/{wildcards.sample}/maker 
		"""


rule collect_gff_annotations:
	input:
		index_path = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.maker.output/{{sample}}_master_datastore_index.log'
	output:
		merged_gff = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.gff'
		
	
	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin && \
		gff3_merge -s -d {input.index_path} > \
			{output.merged_gff}
		"""

rule get_annotation_fasta:
	input:
		index_path = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.maker.output/{{sample}}_master_datastore_index.log'
	output:
		proteins = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.maker.proteins.fasta',
		transcripts = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.maker.transcripts.fasta'

	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin && \
		fasta_merge -d {input.index_path} && \
			mv ./{wildcards.sample}.all.maker.proteins.fasta {output.proteins} && \
			mv ./{wildcards.sample}.all.maker.transcripts.fasta {output.transcripts}
		"""

rule rename_annotation_names:
	input:
		merged_gff = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.gff',

	output:
		map_ids = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.map'
	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin && \
		maker_map_ids --prefix {wildcards.sample}- --justify 6 {input.merged_gff} > {output.map_ids}
		"""

rule rename_proteins_transcripts:
	input:
		proteins = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.maker.proteins.fasta',
		transcripts = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.maker.transcripts.fasta',
		merged_gff = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.gff',
		map_ids = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.map'
	output:
		proteins = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.maker.renamed.proteins.fasta',
		transcripts = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.maker.renamed.transcripts.fasta',
		merged_gff = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.renamed.gff',
	conda:
		"maker"
	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin && \
		cp {input.proteins} ./proteins.fa && \
		cp {input.transcripts} ./transcripts.fa && \
		cp {input.merged_gff} ./annotations.gff && \
		map_gff_ids {input.map_ids} ./annotations.gff && \
		map_fasta_ids {input.map_ids} ./proteins.fa && \
		map_fasta_ids {input.map_ids} ./transcripts.fa && \
		mv ./annotations.gff {output.merged_gff} && \
		mv ./transcripts.fa {output.transcripts} && \
		mv ./proteins.fa {output.proteins}
		"""

rule filter_annotations_default:
	input:
		merged_gff = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.renamed.gff'
	output:
		merged_gff = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.renamed.default.gff'
	conda:
		"maker"
	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin && \
		quality_filter.pl -d {input.merged_gff} > {output.merged_gff}
		"""


rule blast_uniprot:
	input:
		proteins = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.maker.renamed.proteins.fasta',
	output:
		f'{base_dir}/annotation/{{sample}}/functional/uniprot/{{sample}}.all.maker.renamed.maker2uniprot.blastp',
	log:
		'logs/{sample}/annotation/uniprot_blastp.log'
	conda:
		"maker"
	threads:
		40
	shell:
		"""
		blastp -db {config[uniprot_db]} \
		-query {input.proteins} \
		-out {output} \
		-evalue 0.000001 \
		-outfmt 6 \
		-max_hsps 1 \
		-num_alignments 1 \
		-seg yes \
		-lcase_masking \
		-soft_masking true \
		-num_threads {threads} 2> {log}
		"""

rule annotate_functional_uniprot:
	input:
		merged_gff = f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.renamed.default.gff',
		uniprot_blastp = f'{base_dir}/annotation/{{sample}}/functional/uniprot/{{sample}}.all.maker.renamed.maker2uniprot.blastp',
	output:
		f'{base_dir}/annotation/{{sample}}/maker/{{sample}}.all.renamed.default.uniprot.gff'
	conda:
		"maker"
	shadow:
		"shallow"
	shell:
		"""
		export PATH=$PATH:/home/scruz/software/maker/bin && \
		maker_functional_gff {config[uniprot_db]} {input.uniprot_blastp} {input.merged_gff} > {output}
		"""

