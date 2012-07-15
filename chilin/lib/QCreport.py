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
from jinja2 import Environment, FileSystemLoader

#-----------------------
from qc import RawQC,MappingQC,PeakcallingQC,AnnotationQC


# constants
# ------------------------------------

def template_parser():
	tag = {}
	tag['prefix_datasetid'] = 'dataset1277'
	tag['fastqc_check'] = 1
	tag['mapping_check'] = ''
	file = './template/template.tex'
	env = Environment(
			loader=FileSystemLoader('/Users/Samleo/mybin/chilin/chilin/lib/'),
			block_start_string = '\BLOCK{',
			block_end_string = '}',
			variable_start_string = '\VAR{',
			variable_end_string = '}',
			comment_start_string = '\#{',
			comment_end_string = '}',
			line_statement_prefix = '%-',
			line_comment_prefix = '%#',
			trim_blocks = True,
			autoescape = False,
			)
	temp = env.get_template(file)
	print temp.render(prefix_datasetid = tag['prefix_datasetid'],fastqc_check = tag['fastqc_check'],mapping_check = tag['mapping_check'])


def main():
	username = 'meisl'
	outputdir = '/mnt/Storage/home/meisl/mybin/chilin/QCreport/lib/testdata/'
	inputpath = '/mnt/Storage/home/meisl/mybin/chilin/QCreport/lib/QCresult/'
	datasetid = '1277'
	RawQC('sdfsdf').run('yulinyang')
	MappingQC().run()
	PeakcallingQC().run()
	AnnotationQC().run()
#	template_parser()
	
if __name__ == '__main__':
	main()
		
	
	
