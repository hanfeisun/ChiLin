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

#------------------------
from QCreport.Basic_qc import *
from QCreport.fastqc import *
from QCreport.Mapping_qc import *
from QCreport.PeakCalling_qc import *
from QCreport.Ceas_qc import *
from QCreport.Replicate_qc import *
from QCreport.Conservtion_qc import *
from QCreport.motif_qc import *


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


def read_config():
	print 'read_config(): read configfile'



def check_conf_validity():
	print 'check_conf_validity(): check configfile validation'
	
def create_conf():
	print '''create_conf(): Create configfile. The function will search zoho sheet to add some basic information.
	 '''
	
def QC_report():
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
	conservation_info()
	motif_info()
	
	
def main():
	read_config()
	check_conf_validity()
	create_conf()
	QC_report()
	
if __name__ == '__main__':
	main()
		
	
	