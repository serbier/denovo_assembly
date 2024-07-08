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

def get_read_list_per_sample(wildcards):
	read_file = sample_units[sample_units['sample_name'] == wildcards.sample]['read_file_path'].tolist()
	print(read_file)
	return read_file

def get_corrected_reads(wildcards):
	#return "/catalogue/InterspecificCommonBeanDenovoGenomeAssembly/1.Data/corrected/{sample}/correction/NECAT/1-consensus/cns_final.fasta.gz".format(sample = wildcards.sample)
	return "/home/scruz/workflows/ont-basecalling-dna/correction/{sample}/{sample}_corrected.fastq.gz".format(sample = wildcards.sample)



def get_read_depth_cutoffs(wildcards):
	depths = pd.read_csv(config['read_depth_path'], sep='\t')
	depths.set_index('sample', drop=False, inplace=True)

	data = depths.loc[wildcards.sample]
	cutoffs = {
		'l' :data.l,
		'm': data.m,
		'h': data.h	}
	return cutoffs

wildcard_constraints:
    sample = "|".join(sample_units.index),
