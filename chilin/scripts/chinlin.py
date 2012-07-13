#!/usr/bin/python
"""
Last modified
# Time-stamp: <2012-07-11 11:59:42 qinq >
# author: qinq
# status: test

"""
#-------------module--------------
from subprocess import Popen, call
from optparse import OptionParser
import ConfigParser
import string
import re
import logging
import sys
import glob
from os.path import join, exists, getmtime
from jinja2 import Environment, FileSystemLoader, TemplateNotFound
import QC

#----------collective class-------
class FindPath(Check):
    """Find Path according to user filled Meta
    information Table"""
    
    def FindTemporary(self):
        """
        get temporary file in the pipeline
        """
	print "write temporary file name and return to\
               DC or QC"
	print "Write Output file name and Desired Shell Output"
        print "Write Shell Output into temporary file and return path"
    def FindOutput(self):
        print "Simply return output path"



class TemplateParser(Check, FindPath):
    """Load in the QC and DC summary
    template and Write Output in"""
    def Loader(self):
        print "Load in filesystem"
        try:
            env = Environment(loader=FileSystemLoader("chilin/template"))
            template = env.get_template('DA.txt')
            print template.render()
        except TemplateNotFound:
            print "No template folder"

    def DC(self):
        print "Load DC.template and Write in DC output"

    def QC(self):
        print "Load QC template and write into DC output"

class Check(Config, Option, Meta):
    """read in Options, Meta and Config file 
    """
    def __init__(self, conf):
        self._conf = conf

    def ReadConf(self):
        """
        Read configuration and parse it into Dictionary
        """
        print "read and parse the conf"

    def CheckConf(self):
        """
        Check the Meta configuration 
        if up to our definition
        """
        pass

    def DependencyCheck(self):
        """
        check dependency
        """
        print "Check all the dependency program in PATH"
        print "check the input file path"
        print "check output path"
        print "check the external program path"
        print "name the output result and temporary file name"
        print "check shell command line"
        return # True or false for DcController
    def NameOutput(self):
        print "According to the config and set the universal outputname"

class Log(FindPath):
    def __init__(self, option):
        """
        Universal log format
        time, execute, Shell, DC and QC string
        plus time consumed"
        write in the desired log format
        """
        self.log = logging

    def _Timer(self):
        print "Calculate time for each step"
        
    def info(exp):
        print "write running log"
        print "write in shell command"
        print "write accident during pipeline"
        print "_Timer for each step"
        print "Write the path for temporary and output result"

class Bowtie(Check, Log, FindPath):
    """Bowtie DC and QC step"""
    
    def __init__(self):
        print "Get Config and option information Check Class"
        print "Write in Log"
        
    def _format(self):
        print "Get sra or other format into Bowtie Input"
        
    def _Run(self):
        print "Run Bowtie in the Option set command"
        return # true or false for DcController to continue or stop

    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"

class MACS(Check, Log, FindPath, Bowtie):
    """ MACS step, separately and merge for sorted bam
    for peaks calling"""
    
    def __init__(self):
        print "Get MACS related options from options"
        self.options = Options.MACS
        
    def _format(self):
        print "Use samtools Convert the sam to sorted bam"
	print "Convert to tdf format by IGV tools"

    def _Run(self):
        print "call External Macs program"
	print "Print Log into the universal Log"
        
    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"

class CEAS(Check, MACS, Log, FindPath):
    def __init__(self):
        """Get CEAS dependency info from
        Check Class"""
        print self.option
    def _format(self):
        print "Get the random peaks number for speeding the CEAS"

    def _Run(self):
        print "Run the Dependency Program for return FindPath string"
        
    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"

class Seqpos(Check, MACS, Log, FindPath):
    def __init__(self):
	print self.option
    
    def _format(self):
        print "Get the top peaks number for motif analysis"
        
    def _Run(self):
        print "Run the default setting program"
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"

class Replicates(Check, Log, MACS, FindPath):
    def __init__(self):
	print "if replicates, Do this Step"

    def _format(self):
        print "use bedtools to get the desired input"
        print "use bedGraphToBigwiggle to generate Bigwiggle"
    def _Run(self):
        print "Run the venn_diagram"
        print "Run the correlation diagram"
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"

class Conserv(Options, Log, FindPath):
    def _format(self):
	print "Set input Peaks number for conservation Plot"
    def _run(self):
        print "Run the External conservation plot"
    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"

    pass

class GUI(DcController):
    pass

class CistromeAPI(DcController):
    pass

class DcController(Config, Check):
    def _StepControl(self):
        print "Set logic dictionary value for Step"

    def main(self):
        """for running each step of pipeline
        """
        print "run each step in control"
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"
        print "..."

##Options
        # Index-genome: Bowtie
        # m : Bowtie mapping limit
        # model or not : MACS shift-size
        # CEAS : Peaks number
        # Refgene
        # Conservation: Peaks number
        #
        
if __name__ == "__main__":
    main()
