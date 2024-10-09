rule quast:
	input:
		decontaminated_draft = f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.decontaminated.fa",
	output:
		directory(f"{base_dir}/quality/quast/{{sample}}")
	threads:
		20
	conda:
		"quast"
	shell:
		"""
		quast.py -o {output} \
		-t {threads} \
		-r {config[quast_ref]} \
		{input.decontaminated_draft} 
		"""

rule run_busco:
	input:
		f"{base_dir}/quality/decontamination/{{sample}}/{{sample}}.decontaminated.fa",
	output:
		short_json=f"{base_dir}/quality/BUSCO/{{sample}}/short_summary.json",
		short_txt=f"{base_dir}/quality/BUSCO/{{sample}}/short_summary.txt",
		full_table=f"{base_dir}/quality/BUSCO/{{sample}}/full_table.tsv",
		miss_list=f"{base_dir}/quality/BUSCO/{{sample}}/busco_missing.tsv",
		dataset_dir=temp(directory(f"{base_dir}/quality/BUSCO/{{sample}}/busco_downloads"))
	log:
		'logs/{sample}/BUSCO/busco.log'
	params:
		mode="genome",
		lineage="fabales_odb10",
		# optional parameters
		extra="",
	threads:
		40
	wrapper:
		"file:///home/scruz/software/snakemake-wrappers/bio/busco"
