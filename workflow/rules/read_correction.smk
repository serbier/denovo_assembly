rule get_read_list:
	input:
		reads = get_read_list_per_sample,
	output:
		read_list_file = "results/{sample}/correction/read_list.txt"
	run:
		with open(output.read_list_file, "w") as write_file:
			for line in input.reads:
				print(line, file=write_file)

rule produce_necat_config:
	input:
		read_list = "results/{sample}/correction/read_list.txt",
	output:
		"results/{sample}/correction/config.txt"
	params:
		genome_size = config['organism']['genome_size'],
		min_read_l = config['correction']['min_read_l'],

	threads:
		resources['NECAT']['threads']
	shell:
		"""
			python workflow/scripts/get_necat_config.py \
				get_config {wildcards.sample} \
				{params.genome_size} \
				{input.read_list} \
				{params.min_read_l} \
				{threads} \
				{output}
		"""

rule run_correction:
	input:
		config = "results/{sample}/correction/config.txt"
	output:
		directory("results/{sample}/correction/NECAT")
	threads:
		resources['NECAT']['threads']
	log:
		"results/{sample}/correction/correction.log"
	shell:
		"""
			{config[correction][NECAT_path]} correct \
			{input.config} 1> {log} &&
			mv {wildcards.sample} {output}
		"""