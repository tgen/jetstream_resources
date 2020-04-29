#!/usr/bin/env python3

import argparse
import logging as log
import shutil
from sys import exit
from sys import argv
import GTF


class UniqueStore(argparse.Action):
	"""
	class grabbed from stackOverflow (2018-11-08)
	https://stackoverflow.com/questions/23032514/argparse-disable-same-argument-occurences
	Thanks To the Community
	We override the function __call__ from argparse to check if one option is given more than once
	"""

	def __call__(self, parser, namespace, values, option_string):
		if getattr(namespace, self.dest, self.default) is not None:
			parser.error(option_string + " appears several times.  Please modify your options.")
		setattr(namespace, self.dest, values)

def check_if_executable_in_path(list_executables):
	for executable in list_executables:
		if shutil.which(executable) is None:
			exit(str(executable) + "  NOT IN PATH ; Aborting;")


def make_parser_args():
	parser = argparse.ArgumentParser()
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')

	required.add_argument('--gtf',
	                      required=True,
	                      action=UniqueStore,
	                      help='full or relative path to an __Uncompressed__ GTF')
	required.add_argument('--out',
	                      required=True,
	                      action=UniqueStore,
	                      help='name of the output file (full or relative path)')

	return parser


def main(args):
	p = GTF.dataframe(args['gtf'])  ## GFT.dataframe returns a pandas.core.data.DataFrame
	with open(args['out'], 'w') as wo:
		for i in range(len(p)):
			wo.write("{}\t{}\t{}\t{}___{}\t{}\t{}\n".format(p['seqname'][i], p['start'][i], p['end'][i], p['gene_id'][i], p['gene_name'][i],p['gene_biotype'][i], p['strand'][i]))



if __name__ == '__main__':

	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)
	list_executables = ['python3']  ## list of executables required to be in path
	check_if_executable_in_path(list_executables)

	##TODO check if ALL intended python packages are present in python3

	# capture command line
	cmdline = ' '.join(argv)
	log.info(' '.join(["Command Line captured: ", cmdline]))

	# capture arguments
	parser = make_parser_args()
	args = vars(parser.parse_args())  # vars() function returns a dictionary of key-value arguments
	log.info(str(args))
	
	try:
		main(args)
	except Exception as e:
		log.error("ERROR: Exception got Raised; Check you inputs and/or options\n; {}".format(e))
		exit(1)
	exit()

