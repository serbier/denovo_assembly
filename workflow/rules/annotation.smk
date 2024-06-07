rule TE_annotation_scaffolded:
	input:
		scaffolded_draft = 'results/{sample}/scaffolding/{sample}_scaffolded.fasta'
	output:
		te_annotation = '{basedir}/results/{{sample}}/annotation/TE/{{sample}}_TE_annotation_NGSEP.out'.format(basedir=config['basedir'])
	log:
		'{basedir}/results/{{sample}}/annotation/TE/{{sample}}_TE_annotation_NGSEP.log'.format(basedir=config['basedir'])
	conda:
		"denovo"
	threads:
		30
	shell:
		"""
		java -Xmx50G  -jar {config[NGSEP]} TransposonsFinder \
			-i {input.scaffolded_draft} \
			-o {output.te_annotation} \
			-d {config[TE_DB]} \
			-r 2 \
			-t {threads} 2> {log}
		"""

rule convert_NGSEP_TE_to_GFF:
	input:
		te_annotation = '{basedir}/results/{{sample}}/annotation/TE/{{sample}}_TE_annotation_NGSEP.out'.format(basedir=config['basedir'])
	output:
		te_annotation = '{basedir}/results/{{sample}}/annotation/TE/{{sample}}_TE_annotation_NGSEP.gff'.format(basedir=config['basedir'])
	shell:
		"""
		{config[NGSEP_TE2GFF_script]} {input.te_annotation} > \
			{output.te_annotation}
		"""

rule set_maker_config_files:
	input:
		scaffolded_draft = 'results/{sample}/scaffolding/{sample}_scaffolded.fasta',
		te_annotation = '{basedir}/results/{{sample}}/annotation/TE/{{sample}}_TE_annotation_NGSEP.gff'.format(basedir=config['basedir'])


	output:
		bops = '{basedir}/results/{{sample}}/annotation/maker/maker_bopts.ctl'.format(basedir=config['basedir']),
		evm = '{basedir}/results/{{sample}}/annotation/maker/maker_evm.ctl'.format(basedir=config['basedir']),
		exe = '{basedir}/results/{{sample}}/annotation/maker/maker_exe.ctl'.format(basedir=config['basedir']),
		opts = '{basedir}/results/{{sample}}/annotation/maker/maker_opts.ctl'.format(basedir=config['basedir'])
	shadow:
		'shallow'

	conda:
		"tools2"
	shell:
		"""
		{config[maker]} -CTL && \
		awk -v est_path="{config[transcripts]}" \
			-v protein_path="{config[protein]}" \
			-v genome_path="{input.scaffolded_draft}" \
			-v te_path="{input.te_annotation}" -f {config[opts_maker_script]} \
			./maker_opts.ctl > {output.opts} && \
		mv ./maker_bopts.ctl {output.bops} && \
		mv ./maker_evm.ctl {output.evm} && \
		mv ./maker_exe.ctl {output.exe} 
		"""


rule run_maker:
	input:
		bopts = '{basedir}/results/{{sample}}/annotation/maker/maker_bopts.ctl'.format(basedir=config['basedir']),
		evm = '{basedir}/results/{{sample}}/annotation/maker/maker_evm.ctl'.format(basedir=config['basedir']),
		exe = '{basedir}/results/{{sample}}/annotation/maker/maker_exe.ctl'.format(basedir=config['basedir']),
		opts = '{basedir}/results/{{sample}}/annotation/maker/maker_opts.ctl'.format(basedir=config['basedir'])
	output:
		index = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.maker.output/{{sample}}_scaffolded_master_datastore_index.log'.format(basedir=config['basedir'])

	conda:
		"tools2"
	threads:
		22
	log:
		'{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.maker.log'.format(basedir=config['basedir'])
	shell:
		"""
		{config[mpiexec]} -n {threads} {config[maker]} -c 2 -q {input.opts} \
		 {input.bopts} \
		 {input.exe} 2> {log} && \
		 mv ./{wildcards.sample}_scaffolded.maker.output {config[basedir]}/results/{wildcards.sample}/annotation/maker
		"""


rule collect_gff_annotations:
	input:
		index_path = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.maker.output/{{sample}}_scaffolded_master_datastore_index.log'.format(basedir=config['basedir'])

	output:
		merged_gff = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.all.gff'.format(basedir=config['basedir'])
	
	shell:
		"""
		{config[gff3_merge]} -s -d {input.index_path} > \
			{output.merged_gff}
		"""

rule get_annotation_fasta:
	input:
		index_path = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.maker.output/{{sample}}_scaffolded_master_datastore_index.log'.format(basedir=config['basedir'])
	output:
		proteins = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.all.maker.proteins.fasta'.format(basedir=config['basedir']),
		transcripts = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.all.maker.transcripts.fasta'.format(basedir=config['basedir'])
	shell:
		"""
		{config[fasta_merge]} -d {input.index_path} && \
			mv ./{wildcards.sample}_scaffolded.all.maker.proteins.fasta {output.proteins} && \
			mv ./{wildcards.sample}_scaffolded.all.maker.transcripts.fasta {output.transcripts}
		"""

rule rename_annotation_names:
	input:
		merged_gff = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.all.gff'.format(basedir=config['basedir']),
		proteins = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.all.maker.proteins.fasta'.format(basedir=config['basedir']),
		transcripts = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.all.maker.transcripts.fasta'.format(basedir=config['basedir'])

	output:
		map_ids = '{basedir}/results/{{sample}}/annotation/maker/{{sample}}_scaffolded.all.map'.format(basedir=config['basedir'])
	shell:
		"""
		{config[maker_map_ids]} --prefix {wildcards.sample} --justify 6 {input.merged_gff} > {output.map_ids} && \
		{config[map_gff_ids]} {output.map_ids} {input.merged_gff} && \
		{config[map_fasta_ids]} {output.map_ids} {input.proteins} && \
		{config[map_fasta_ids]} {output.map_ids} {input.transcripts}
		"""

