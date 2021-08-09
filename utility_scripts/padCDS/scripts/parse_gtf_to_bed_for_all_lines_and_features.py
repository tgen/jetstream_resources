#!/usr/bin/env python3

import argparse
import logging
import shutil
from os import system, popen
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


def add_field_to_GTF(field, gtf_file):
	cmd = 'echo -n $( head -n 3 ' + str(gtf_file) + ' | tail -n 1  ) | tail -c 3 | grep -c ";"'
	logger.debug(cmd)
	c = popen(cmd).read()
	logger.debug("count of semi-colon: " + str(c.strip()))
	if c.strip() == "1":
		logger.debug("not adding semi-colon".upper())
		cmd = 'sed -i \'s/$/ '+field+' \".\";/\''
	else:
		logger.debug(" adding semi-colon".upper())
		cmd = 'sed -i \'s/$/; '+field+' \".\";/\''
	cmd = cmd + " " + gtf_file
	logger.debug(str(cmd))
	try:
		system(cmd)
	except Exception as E:
		logger.error(e)
	logger.info("field update done for: "+str(field))


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

	# --------------------------------------------------------------------------
	# READING GTF with GTF.py
	# --------------------------------------------------------------------------
	try:
		gtf_file=args['gtf']
		p = GTF.dataframe(gtf_file)  ## GTF.dataframe returns a pandas.core.data.DataFrame
	except Exception as e:
		logger.error("ERROR: in reading GTF\n; {}".format(e))
		exit(1)

	# --------------------------------------------------------------------------
	# INIT VARIABLES
	# --------------------------------------------------------------------------
	_Expected_Eight_Columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
	_Expected_Other_Columns = [ 'gene_id', 'gene_name', 'gene_version', 'gene_source', 'gene_biotype', 'transcript_id', 'transcript_version', 'transcript_name', 'transcript_source', 'transcript_biotype', 'tag', 'transcript_support_level', 'exon_number', 'exon_id', 'exon_version', 'protein_id', 'protein_version' ]

	logger.debug(p.index)
	logger.debug(p.columns)
	logger.info("-"*50)
	logger.info("Printing the fields (aka columns) needed in GTF to make the BED:")
	logger.info(_Expected_Eight_Columns + _Expected_Other_Columns)
	logger.info("-"*50)

	# --------------------------------------------------------------------------
	# CHECK FIELDS
	# --------------------------------------------------------------------------
	mandatory_columns_count = 0
	mandatory_columns_list = []
	other_fields_list = []
	for col in p.columns:
		if col in _Expected_Eight_Columns:
			mandatory_columns_count +=1
			mandatory_columns_list.append(col)
		other_fields_list.append(col)
	if mandatory_columns_count != len(_Expected_Eight_Columns):
			raise ValueError("MISSING MANDATORY COLUMNS in GTF: {} ".format(';'.join([ str(col) for col in _Expected_Eight_Columns if col not in mandatory_columns_list ])) )
	logger.info("Mandatory Expected Fields Check out OK")
	logger.info("Testing if any other missing field exists or are extra or with different expected names ...")

	unexpected_fields_list = []
	expected_fields_list = []
	for field in other_fields_list:
		if field not in _Expected_Other_Columns and field not in _Expected_Eight_Columns:
			unexpected_fields_list.append(field)
		else:
			expected_fields_list.append(field)

	logger.warning("The Following Fields are NOT added to BED; Check if they should be used and if so, check if they might be named differently: {} ".format(unexpected_fields_list))

	# --------------------------------------------------------------------------
	# UPDATE GTF IF NEEDED
	# --------------------------------------------------------------------------
	read_GTF_again=False
	for field in _Expected_Other_Columns + _Expected_Eight_Columns:
		if field not in expected_fields_list:
			if not read_GTF_again:
				shutil.copy(gtf_file, gtf_file+"upd.gtf")
				gtf_file = gtf_file+"upd.gtf"
			add_field_to_GTF(field, gtf_file)
			read_GTF_again=True

	if read_GTF_again:
		logger.info("processing GTF2BED for << {} >> updated GTF...".format(gtf_file))
		p = GTF.dataframe(gtf_file)

	# --------------------------------------------------------------------------
	# writing to output file
	# --------------------------------------------------------------------------
	with open(args['out'], 'w') as wo:
		logger.info("writing out BED file ... {}".format(args['out']))
		## writing HEADER line to output file
		wo.write("\t".join(
			['##seqname', 'start', 'end', 'gene_id__gene_name', 'score', 'strand', 'frame', 'gene_version', 'gene_source', 'gene_biotype', 'transcript_id', 'source', 'feature', 'transcript_version', 'transcript_name', ' transcript_source', 'transcript_biotype', 'tag',
			 'transcript_support_level', 'exon_number', 'exon_id', 'exon_version', 'protein_id', 'protein_version']) + "\n")
		## writing VALUE lines to output file
		try:
			for i in range(len(p)):
				wo.write("{}\t{}\t{}\t{}___{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(p['seqname'][i], int(p['start'][i]) - 1, p['end'][i], p['gene_id'][i], p['gene_name'][i], p['score'][i], p['strand'][i], p['frame'][i], p['gene_version'][i], p['gene_source'][i], p['gene_biotype'][i], p['transcript_id'][i], p['source'][i], p['feature'][i], p['transcript_version'][i], p['transcript_name'][i], p['transcript_source'][i], p['transcript_biotype'][i], p['tag'][i], p['transcript_support_level'][i], p['exon_number'][i], p['exon_id'][i], p['exon_version'][i], p['protein_id'][i], p['protein_version'][i] )
				)
		except IOError as IOE:
			logger.error("ERROR: in Writing data\n; {}".format(IOE))
			exit(2)
		except ValueError as VE:
			logger.error("ERROR: in Writing data\n; {}".format(VE))
			exit(2)
		except Exception as E:
			logger.error("ERROR: in Writing data\n; {}".format(E))
			exit(1)


if __name__ == '__main__':

	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	logging.basicConfig(format=FORMAT_LOGGING, level=logging.INFO)
	logger = logging.getLogger(__name__)

	list_executables = ['python3']  ## list of executables required to be in path
	check_if_executable_in_path(list_executables)

	# capture command line
	cmdline = ' '.join(argv)
	logger.info(' '.join(["Command Line captured: ", cmdline]))

	# capture arguments
	parser = make_parser_args()
	args = vars(parser.parse_args())  # vars() function returns a dictionary of key-value arguments
	logger.info(str(args))

	try:
		main(args)
	except Exception as e:
		logger.error("ERROR: Exception got Raised; Check you inputs and/or options\n; {}".format(e))
		exit(1)
	exit()
