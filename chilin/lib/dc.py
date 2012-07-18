from subprocess import call
from jinja2 import Environment, PackageLoader, TemplateNotFound
from ConfigParser import ConfigParser, SafeConfigParser
import logging
import string
import sys
import glob
import datetime
import os

jinja_env = Environment(loader = PackageLoader('chilin', 'template'),
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
cf = ConfigParser()
safecf = SafeConfigParser()
class PipePreparation:
    """read in Options set by optparse
    """
    def __init__(self, ChiLinconfPath = ''):
        self.ChiLinconfPath = ChiLinconfPath
        self.ChiLinconfigs = {}
        print "read in options from command line"

    def _readconf(self):
        """
        Read configuration and parse it into Dictionary
        """
        cf.read(self.ChiLinconfPath)
        for sec in cf.sections():
            secName = string.lower(sec)
            for opt in cf.options(sec):
                print opt
                optName = string.lower(opt)
                self.ChiLinconfigs[secName + '.' + optName] = string.strip(cf.get(sec, opt))

    def checkconf(self):
        """
        Check the Meta configuration
        if up to our definition
        """
        self._readconf()

        if not os.path.exists(self.ChiLinconfigs['bowtie.bowtie_main']):
            print "bowtie program dependency has passed"
            return False
        if not os.path.exists(self.ChiLinconfigs['macs.macs_main']):
            print "macs2 program dependency has passed"
            return False
        return True

class PathFinder:
    """prepare path for each step"""
    def __init__(self, NameConf='', rep='', datasetid=''):
        self.cf = ConfigParser()
        self.NameConf = NameConf
    def _readconf(self):
        self.cf.read(NameConf)
    def bowtiefilepath(self):
        print self.cf.sections()
        cmd = 'mkdir %s' % self.bowtiefolder
        call(cmd, shell = True)
        return
    def macs2filepath(self):
        print 'macs2filepath'
        return 'path'
    def venn_corfilepath(self):
        return "path"
    def ceasfilepath(self):
        return "path"
    def conservfilepath(self):
        return 'path'
    def qcfilepath(self):
        return "qcfilepath"

class LogWriter:
    def __init__(self, logfile = 'log'):
        """
        Universal log format
        time, execute, Shell, DC and QC string
        plus time consumed
        """
        self.logfile = logfile

    def record(self, logcontent):
        dt = datetime.datetime.now()
        with open(self.logfile) as f:
            f.write(dt.strftime('%Y-%m-%d-%H:%M:%S') + logcontent)

class PipeController(object):
    def __init__(self, Options = ""):
        """
        # read in Options from command line
        # Get template and conf information
        """
        print "Dc control prepare"

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
        cmd = 'cp {0} {1}'
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

class PipeBowtie(PipeController):
    """Bowtie DC and QC step"""
    
    def __init__(self):
        super(PipeBowtie, self).__init__()
        self.bowtie_main = self.parser.get("bowtie", "bowtie_main")
        self.gene_index = self.parser.get("bowtie", "bowie_genome_index")
        
    def _format(self):
        print "Get sra or other format into Bowtie Input"
        print "Write in Log"
        
    def _run(self):
        cmd = '{0} -S {1} -m {2} {3} {4} {5}'
        print "Run Bowtie in the Option set command"
        print "Write in Log"
        return # true or false for PipeController to continue or stop\
                # and path for Next step

    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Write into the template"

class PipeMACS2(PipeController):
    """ MACS step, separately and merge for sorted bam
    for peaks calling"""
    
    def __init__(self):
        super(PipeMACS2, self).__init__()
        
    def _format(self):
        return 

    def _run(self):
        cmd = '{0} callpeak {1} -B -q 0.01 --nomodel --shift-size {2} ' + \
              '-t {3} -c {4} -n {5}'
               
            
        cmd = cmd.format(self.macs2_main,
                         self.genome_option,
                         self.shiftsize,
                         self.treat_bam,
                         self.control_bam,
                         self.macsname)
        
        return cmd

    def summary(self):
        call(self._run())


class PipeVennCor(PipeController):
    def __init__(self, OptionMethod = "Mean"):
        """Read in the Controller Dictionary
        to decide whether do this step or not"""
        super(PipeVennCor, self).__init__()
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

class PipeCEAS(PipeController):
    def __init__(self):
        """Get CEAS dependency info from
        Check Class"""
        self.format = []
        super(PipeCEAS, self).__init__()

    def _format(self):
        self.format.append('BED or wiggele or both')
        print "Get the random peaks number for speeding the CEAS"

    def _run(self):

        print "Run the Dependency Program for return FindPath string"
        
    def summary(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Write into the template"

class PipeConserv(PipeController):
    def __init__(self):
        super(PipeConserv, self).__init__()

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

class PipeMotif(PipeController):
    def __init__(self):
        super(PipeMotif, self).__init__()
    
    def _format(self):
        print "Get the top peaks number for motif analysis"
        
    def _Run(self):
        print "Run the default setting program"
        print "Call private _Run"
        print "Extract shell output"
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"
