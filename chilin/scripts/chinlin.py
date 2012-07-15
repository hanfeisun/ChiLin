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
from os.path import join, exists, getmtime
from jinja2 import Environment, FileSystemLoader, TemplateNotFound
from collections import defaultdict

import ConfigParser
import string
import re
import logging
import sys
import glob
import time
import datetime
#import QC


#----------collective class-------

class Check(Config, Options, Meta):
    """read in Options set by optparse
    Meta -> ${DatasetID}Meta.xls
    Name rule file -> NameRule.ini

    return a dictionary as follows:
    PathFindDict = {${DatasetID}_bowtie.sam : 
                   [ path/to/file+universalfilename, True] #logic value for control the step
    """
    def __init__(self, conf):
        self.tmpname = ''
        self.outputname = ''
        self.reportname = ''
        self.cormethod = options.method

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
        print "Check the configuration right or not"

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

    def PathFinder(self):
        """
        According to Output&temporary file name
        and user added Meta data path,
        return a defaultDictionary which reflects
        the replicates or not, good quality or poor quality
        thus, control the step
        """
        print "read in the output contemporary Name from config file"
        print "read in the Meta data path from .xls file"
        print "According to the config and set the universal outputname"

class TemplateParser(Check):
    """Load in the QC and DC summary
    template and Write Output in"""
    def __init__(self):
        self.datatype = Check
        super(TemplateParser, self).__init__()

    def Loader(self):
        print "Load in filesystem"
        try:
            env = Environment(loader=FileSystemLoader("chilin/template"))
            template = env.get_template('DA.txt')
            #template = env.get_template('template.tex')
            print template.render()
        except TemplateNotFound:
            print "No template folder"

    # write in DC, QC template separate variables
    def DcRender(self):
        print "Load DC.template and Write in DC output"

    def QcRender(self):
        print "Load QC template and write into DC output"


class Log(logging):
    def __init__(self, option):
        """
        Universal log format
        time, execute, Shell, DC and QC string
        plus time consumed"
        write in the desired log format
        """
        basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    ) # logging
        self.logfhd = open("log","w")
        self.error   = critical # alias for logging.critical
        self.warn    = warning

    def Timer(self):
        print "Calculate time for each step"
        print "CPU time elapsed: %s s" %(time.clock() - cpustart)
        print "Real time elapsed: %s s" %(datetime.datetime.now()-realstart)
        
    def info(exp):
        self.info(a) # logging.info
        self.logfhd.write(a + '\n')
        self.logfhd.flush()
        if warn:
            self.warn()
        elif error:
            self.error()
        print "write running log"
        print "write in shell command"
        print "write accident during pipeline"
        print "_Timer for each step"
        print "Write the path for temporary and output result"

#---------------DC class---------------

class Bowtie(Check, Log):
    """Bowtie DC and QC step"""
    
    def __init__(self):
        print "Get Config and option information Check Class"
        super(Bowtie, self).__init__()
        
    def _format(self):
        print "Get sra or other format into Bowtie Input"
        print "Write in Log"
        
    def _run(self):
        print "Run Bowtie in the Option set command"
        print "Write in Log"
        return # true or false for DcController to continue or stop

    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"

class MACS(Check, Log):
    """ MACS step, separately and merge for sorted bam
    for peaks calling"""
    
    def __init__(self):
        print "Get MACS related options from options"
        
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
        """Read in the Controller Dictionary
        to decide whether do this step or not"""
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

class Conserv(Check, Log):

    def __init__(self):
        super(Conserv, self).__init__()

    def _format(self):
        print "Set input Peaks number for conservation Plot"
    def _run(self):
        print "Run the External conservation plot"
    def summary(self):
        print "Call private _Run"
        print "Extract shell output" # generate temporary file
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"
        return # information for passing to TemplateParser

class GUI(DcController):
    pass

class CistromeAPI(DcController):
    def __init__(self):
        super(DcController)
    pass

class DcController(Config, Check):

    def _StepControl(self):
        """ Supplement Dictionary logic value 
        for class Check"""
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
