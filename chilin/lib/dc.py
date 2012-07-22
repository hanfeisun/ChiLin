from subprocess import call
from jinja2 import Environment, PackageLoader, TemplateNotFound
from ConfigParser import ConfigParser, SafeConfigParser
import logging
import string
import sys
import glob
import datetime
import os
import chilin

#alias
j = os.path.join
e = os.path.exists

#class for configs
class PipePreparation:
    """read in Options set by optparse
    """
    def __init__(self, ChiLinconfPath = ''):
        self.ChiLinconfPath = ChiLinconfPath
        self.cf = ConfigParser()
        self.ChiLinconfigs = {}

    def readconf(self):
        """
        Read configuration and parse it into Dictionary
        """
        self.ChiLinconfPath
        self.cf.read(self.ChiLinconfPath)
        for sec in self.cf.sections():
            temp = {}
            for opt in self.cf.options(sec):
                optName = string.lower(opt)
                temp[string.strip(optName)] = string.strip(self.cf.get(sec, opt))
            self.ChiLinconfigs[string.lower(sec)] = temp
    def checkconf(self):
        """
        Check the Meta configuration
        if up to our definition
        """
        self.readconf()
        if not e(self.ChiLinconfigs['qc']['fastqc_main']):
            print 'fastqc not exists'
            return False
        if not e(self.ChiLinconfigs['bowtie']['bowtie_main']):
            print "bowtie program dependency has failed"
            return False
        if not e(self.ChiLinconfigs['macs']['macs_main']):
            print "macs2 program dependency has failed"
            return False
        return True

class PathFinder:
    """prepare path for each step"""
    def __init__(self, outputd = '',  datasetid='', treatpath = '', controlpath = '',\
            NameConfPath = os.path.split(chilin.__file__)[0] + '/' + 'db/NameRule.conf'):
        self.cf = SafeConfigParser()
        self.NameConfPath = NameConfPath
        self.Nameconfigs = {}
        self.treat_path = treatpath.split(',')
        self.control_path = controlpath.split(',')
        self.datasetid = datasetid
        self.outputd = outputd

    def _readconf(self):
        '''read in conf and write datasetid information'''
        self.cf.read(self.NameConfPath)
        for sec in self.cf.sections():
            temp = {}
            for opt in self.cf.options(sec):
                optName = string.lower(opt)
                temp[string.strip(optName)] = string.strip(self.cf.get(sec, opt).replace('${DatasetID}', self.datasetid))
                self.Nameconfigs[string.lower(sec)] = temp

    def parseconfrep(self):
        '''only name means not plus the output directory, for legend only
        write in NameRule rep '''
        self._readconf()

        for session in self.Nameconfigs:
            for option in self.Nameconfigs[session]:
                temp = []
                if self.control_path[0] != '':
                    if '${control_rep}' in self.Nameconfigs[session][option]:
                        for control_rep in range(1, len(self.control_path) + 1):
                            if not self.outputd.endswith('/'):
                                temp.append(self.outputd + '/' + self.Nameconfigs[session][option].replace('${control_rep}', str(control_rep)))
                            else:
                                temp.append(self.outputd + self.Nameconfigs[session][option].replace('${control_rep}', str(control_rep)))
                        self.Nameconfigs[session][option] = temp

                if self.treat_path[0] != '':
                    if '${treat_rep}' in self.Nameconfigs[session][option]:
                        for treat_rep in range(1, len(self.treat_path) + 1):
                            if not self.outputd.endswith('/'):
                                temp.append(self.outputd + '/' + self.Nameconfigs[session][option].replace('${treat_rep}', str(treat_rep)))
                            else:
                                temp.append(self.outputd + self.Nameconfigs[session][option].replace('${treat_rep}', str(treat_rep)))

                        self.Nameconfigs[session][option] = temp

# log 
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
        with open(self.logfile, 'wb') as f:
            f.write(dt.strftime('%Y-%m-%d-%H:%M:%S') + logcontent)
            logging.error(logcontent)
            logging.warning(logcontent)

# parent class of DC part
class PipeController(object):
    def __init__(self):
        """
        read in Options from command line
        Get template and conf information
        """
        self.has_run = True
        self.cmd = ''

    def run(self, pipestep = ''):
        """for running each step of pipeline
        """
        shelljudge = call(self.cmd, shell = True)
        print self.chilinconfigs
        if shelljudge == 0:
            self.has_run = True
        else:
            self.has_run = False

    def partition(self):
        self.cmd = 'mv {0} {1}'
        mv = call(cmd, shell = True)
        if mv:
            return True

    def render(self, template = ''):
        """
        write into the DA.txt template
        """
        if self.has_run:
            print "Get variable"
            print "Write into Template"
        else:
            print "Write into log"

# children class
class PipeBowtie(PipeController):
    """Bowtie DC and QC step"""
    def __init__(self, chilinconfigs, nameconfigs):
        super(PipeBowtie, self).__init__()
        self.chilinconfigs = chilinconfigs
        self.nameconfigs = nameconfigs

    def _format(self):
        if 'format' == 'fastq':
            print "Get sra or other format into Bowtie Input"
            print "Write in Log"
        else:
            pass

    def process(self):
        self.run()
        self.cmd  = '{0} -S {1} -m {2} {3} {4} {5}'
        print 'test', self.nameconfigs
        print 'test2', self.chilinconfigs

        self.cmd.format(self.chilinconfigs['bowtie']['bowtie_main'],
                        'sam',
                        max_alignment,
                        index_path,
                        treat_path,
                        outputname)
        call(self.cmd)
        self._format()

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
