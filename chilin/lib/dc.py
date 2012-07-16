from subprocess import Popen, call
from os.path import join, exists
from jinja2 import Environment, FileSystemLoader, TemplateNotFound
from collections import defaultdict
from ConfigParser import ConfigParser

import string
import re
import logging
import sys
import glob
import time
import datetime
import os

class Check(object):
    """read in Options set by optparse
    Meta -> ${DatasetID}Meta.xls
    Name rule file -> NameRule.ini

    return a dictionary as follows:
    PathFindDict = {${DatasetID}_bowtie.sam : 
                   [ path/to/file+universalfilename, True] #logic value for control the step
    """
    def __init__(self, config="", NameRule="", Meta=""):
        self.conf = os.path.join(os.path.split(os.getcwd())[0], 'lib/db/chinlin.ini')
        self.Name = NameRule
        self.Meta = Meta # from options
        print "read in options from command line"

    def ReadConf(self):
        """
        Read configuration and parse it into Dictionary
        """

        print self.conf
        cf = ConfigParser()
        cf.read(self.conf)
        secs = cf.sections()
        print secs
        print cf.options('bowtie')

    def MetaParse(self, MetaPath = ""):
        print MetaPath

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
        print "read in NameRule.ini to initialize output& temporary file name"
        print "Create corresponding folder"
        print "Set path for each file according to Meta.xls output column"
        print "write out step control information using a Dictionary"

class TemplateParser(object):
    """Load in the QC and DC summary
    template and Write Output in"""
    def __init__(self):
        self.datatype = Check
        super(TemplateParser, self).__init__()

    def Loader(self):
        print "Load in filesystem"
        print "Create a temporary file to contain writen summary"
        try:
            env = Environment(loader=FileSystemLoader("chilin/template"))
            template = env.get_template('DA.txt')
            #template = env.get_template('template.tex')
            print template.render()
        except TemplateNotFound:
            print "No template folder"

class Log(object):
    def __init__(self):
        """
        Universal log format
        time, execute, Shell, DC and QC string
        plus time consumed"
        write in the desired log format
        """
        print "Create a Log file"
        print "Customized log content, including shell output,\
                warning and error"

    def Timer(self):
        print "Calculate time for each step"
        
    def info(exp):
        print "write config and meta information in log file"
        print "write shell command in log"
        print "if warning, write in warning"
        print "if error, break, write in error"
        print "if execute command succeed, call Timer to write in running time"

class DcController(object):
    def __init__(self, Options = ""):
        """
        # read in Options from command line
        # Get template and conf information
        """
        print "Dc control prepare"
        self.error = False
        print "DC control start"



    def _StepControl(self):
        """ Supplement Dictionary logic value 
        for class Check"""
        print "Set logic dictionary value from Options and Check class"

    def run(self):
        """for running each step of pipeline
        """
        print "run each step in control"
        print "Call private class bowtie etc. _run"
        print "Extract shell output for render and log"
        if self.error: 
            print "write log and exit"
            return False # logic
        else:
            print "write in warning"
            return True

    def partition(self):
        print "Create Folder up to output&temporary folder"
        print "assign the output&temporary file according to\
               output& temporary folder"

    def render(self, template = ''):
        """
        write into the DA.txt template 
        separately 
        """
        if self.run:
            print "Get variable"
            print "Write into Template"
        else:
            print "Write into log"


class Bowtie(DcController):
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
        return # true or false for DcController to continue or stop\
                # and path for Next step

    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Write into the template"

class MACS(DcController):
    """ MACS step, separately and merge for sorted bam
    for peaks calling"""
    
    def __init__(self):
        self.Model = True
        print "Get MACS related options from options"
        super(Replicates, self).__init__()
        
    def _format(self):
        print "Use samtools Convert the sam to sorted bam"
        print "Convert to tdf format by IGV tools"

    def _run(self):
        if self.Model:
            print "Build Model"
        else:
            print "--shift-size --nomode options on"

        print "call External Macs program"
        print "Print Log into the universal Log"
        
    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Write into the template"

class Replicates(DcController):
    def __init__(self, OptionMethod = "Mean"):
        """Read in the Controller Dictionary
        to decide whether do this step or not"""
        super(Replicates, self).__init__()
        print "if replicates, Do this Step"

    def _format(self):
        print "use bedtools to get the desired input"
        print "use bedGraphToBigwiggle to generate Bigwiggle"
    def _run(self):
        print "Run the venn_diagram"
        print "Run the correlation diagram"
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"
        return "Path for QC"

    def summary(self):
        print "render to template"

class CEAS(DcController):
    def __init__(self):
        """Get CEAS dependency info from
        Check Class"""
        self.format = []
        super(Replicates, self).__init__()

    def _format(self):
        self.format.append('BED or wiggele or both')
        print "Get the random peaks number for speeding the CEAS"

    def _run(self):

        print "Run the Dependency Program for return FindPath string"
        
    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Write into the template"

class Seqpos(DcController):
    def __init__(self):
        super(Replicates, self).__init__()
    
    def _format(self):
        print "Get the top peaks number for motif analysis"
        
    def _Run(self):
        print "Run the default setting program"
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"


class Conserv(DcController):
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



# TO-DO
class GUI(DcController):
    pass

class CistromeAPI(DcController):
    pass
