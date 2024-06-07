import os, sys, re, argparse


class main():
	def __init__(self):
		parser = argparse.ArgumentParser(
		description='Get NECAT configuration file',
		usage='''get_necat_config.py <command> [<args>]
		Available commands are:
		get_config    Get configuration file according to parameters

		''')

		parser.add_argument('command', help='Subcommand to run')
		args = parser.parse_args(sys.argv[1:2])

		if not hasattr(self, args.command):
			print('Unrecognized command')
			parser.print_help(file=sys.stderr)
			exit(1)
		getattr(self, args.command)(sys.argv[2:])

	def get_config(self, argv):
		parser = argparse.ArgumentParser(description="Get configuration file according to parameters")
		parser.add_argument("project_name", help="Name of the project")
		parser.add_argument("genome_size", help="Genome size of organism in bp")
		parser.add_argument("read_list_path", help="Path of read list file")
		parser.add_argument("min_read_l", help="Minimum read length to correct")
		parser.add_argument("threads", help="Thread number to use")
		parser.add_argument("out_file", help="Path where config file will be stored")

		args = parser.parse_args(argv)

		config_txt = """PROJECT={project_name}
ONT_READ_LIST={read_list_path}
GENOME_SIZE={genome_size}
THREADS={threads}
MIN_READ_LENGTH={min_read_l}
PREP_OUTPUT_COVERAGE=40
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true""".format(project_name =args.project_name,
	read_list_path =args.read_list_path,
	genome_size =args.genome_size,
	min_read_l=args.min_read_l,
	threads =args.threads)

		with open(args.out_file,"w+") as f:
			f.writelines(config_txt)

if __name__ == "__main__":
	main()