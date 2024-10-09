import pandas as pd
from yaml import safe_load


# Snakemake configuration
configfile: "config/config.yaml"

with open(config["resources_config"], "r") as f:
    resources = safe_load(f)

config['organism']['genome_size_mb'] = int(config['organism']['genome_size']/1e6)


sample_units = pd.read_table(config["sample_units"], sep='\t')
sample_units.set_index('sample_name', drop=False, inplace=True)

base_dir = '/opt/bean_fast/ONT'

""" Input function used in raw_reads_processing/merge_raw_reads rule """

def get_read_list_per_sample(wildcards):
	read_file = sample_units[sample_units['sample_name'] == wildcards.sample]['read_file_path'].tolist()
	return read_file

def get_raw_reads(wildcards):
	reads = get_read_list_per_sample(wildcards)
	if len(reads) > 1:
		return f"{base_dir}/basecalling/{wildcards.sample}/merged.fastq.gz"
	else:
		return reads[0]
		
def get_read_depth_cutoffs(wildcards):
	depths = pd.read_csv(config['read_depth_path'], sep='\t')
	depths.set_index('sample', drop=False, inplace=True)

	data = depths.loc[wildcards.sample]
	cutoffs = {
		'l' :data.l,
		'm': data.m,
		'h': data.h	}
	return cutoffs


""" Bring assembly to local disk to avoid
I/O bottleneck """
rule bring_assembly:
	input:
		draft = f'{base_dir}/assembly/flye/{{sample}}/assembly.fasta',
	output:
		draft = temp("data/assembly/{sample}.fasta"),
		index = temp("data/assembly/{sample}.fasta.fai")
	conda:
		"ngs"
	shell:
		"""
		cp {input.draft} {output.draft} && \
		samtools faidx {output.draft}
		"""

wildcard_constraints:
    sample = "|".join(sample_units.index),
