#!/usr/bin/python


# python modules
import os
import sys
import re
from optparse import OptionParser
import ConfigParser
import string
import logging
from subprocess import call

#-----------------------
from .basic_qc import Basic_info,fastqc_info
from .mapping_qc import Mapping_info
from .peak_calling_qc import PeakSummary_info,DHS_info,Velcro_info,Replicate_info
from .function_qc import Motif_info,Ceas_info,Conservation_info


# constants
# ------------------------------------

logfhd = open("log","w")

logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

error   = logging.critical        # function alias
warn    = logging.warning

def info(a):
    logging.info(a)
    logfhd.write(a+"\n")
    logfhd.flush()


#
def datafind(d,pattern):
	print 'The function can be used to find particular file. '


def read_config(configFile):
	print 'read_config(): read configfile'
	"""Read configuration file and parse it into a dictionary.
	In the dictionary, the key is the section name plus option name like: data.data.treatment_seq_file_path.
    """
	configs = {}
	configParser = ConfigParser.ConfigParser()
	if len(configParser.read(configFile)) == 0:
		raise IOError("%s not found!" % configFile)
	for sec in configParser.sections():
		secName = string.lower(sec)
		for opt in configParser.options(sec):
			optName = string.lower(opt)
			configs[secName + "." + optName] = string.strip(configParser.get(sec, opt))
	
	return configs

def laytex_parser(f):
	infile = open(f)
	content = infile.readlines()
	contentr=''
	for i in range(len(content)):
		contentr +=content[i]
	temporate = {}
	QC=['prefix_tex','basic_tex','fastqc_tex','mapping_qc_tex','peak_summary_qc_tex','dhs_qc_tex','velcro_qc_tex','replicate_qc_tex','ceas_qc_tex','conservation_qc_tex','motif_qc_tex']
	for i in range(len(QC)):
		t='<%s>[\w^\W]<%s>' %(QC[i],QC[i])
		rule = re.compile(r'<%s>[\w^\W]*<%s>' %(QC[i],QC[i]))
		lines = rule.findall(contentr)
		if lines:
			cut=len(QC[i])+2
			temporate[QC[i]] = lines[0][cut:-cut]
		else:
			print t
	return temporate


def check_conf_validity():
	print 'check_conf_validity(): check configfile validation'
	
def create_conf(username,datasetid,inputpath,outputdir,configs):
	print 'create_conf(): Create configfile. The function will search zoho sheet to add some basic information.'
	""" Create a config file and output directory"""
	cmd = 'mkdir %s' %outputdir
	call(cmd,shell=True)
	if os.path.exists(outputdir):
		info('output folder bulit')
	else:
		error("can't bulit output folder")
	zoho = open(configs["path.zoho_sheet"],'rU')
	zoho_infile = zoho.readlines()
	for i in range(len(zoho_infile)):
		temp = zoho_infile[i].split("\t")
		if str(temp[0]) == str(datasetid):
			configs["basic.species"] = temp[8].lower()
			cell_type = temp[9].lower()
			cell_line = temp[11].lower()
			factor = temp[16]
			format = temp[25].strip().split(".")
			format = format[-1]
			replicate = len(temp[25].split(","))
#			configs["basic.treatment"] = temp[24]
#			configs["basic.control"] = temp[26]
			print factor
			print replicate
			print format
			break
	
	configs["basic.username"] = username
	configs["basic.dataset_id"] = datasetid
	if inputpath.endswith('/'):
		configs["basic.dataset_path"] = inputpath
	else:
		configs["basic.dataset_path"] = inputpath + '/'
	if outputdir.endswith('/'):
		configs["basic.output_path"] = outputdir
	else:
		configs["basic.output_path"] = outputdir + '/'
	configs["basic.raw_data_format"] = format
	configs["basic.replicate"] = replicate
	configs["basic.cell_type"] = cell_type
	configs["basic.cell_line"] = cell_line
	configs["basic.factor"] = factor
	
	configs["path.ceas_code"] = ''
	configs["path.conservation_code"] = ''
	configs["path.conservation_pdf"] = ''
	configs["path.summary_file"] = ''
		
	return configs


	
def QC_report():
	#-------mapping QC---------
	Basic_info()
	fastqc_info()
	Mapping_info()
	
	#--------Peak calling measurement
	PeakSummary_info()
	DHS_info()
	Velcro_info()
	Replicate_info()
	
	#---------Functional Genomic QC measurement----------
	Ceas_info()
	Conservation_info()
	Motif_info()
	
	
def main():
	configs = read_config("./db/config.conf")
	username = 'meisl'
	outputdir = '/mnt/Storage/home/meisl/mybin/chilin/QCreport/lib/testdata/'
	inputpath = '/mnt/Storage/home/meisl/mybin/chilin/QCreport/lib/QCresult/'
	datasetid = '1277'
	check_conf_validity()
	configs = create_conf(username,datasetid,inputpath,outputdir,configs)
	print configs
	QC_report()
	latex_tex=laytex_parser('./template/latex_template.tex')
	print latex_tex
	
if __name__ == '__main__':
	main()
		
	
	
