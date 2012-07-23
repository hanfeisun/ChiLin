from subprocess import call, check_output
from jinja2 import Environment, PackageLoader, TemplateNotFound
from ConfigParser import ConfigParser, SafeConfigParser
import logging
import string
import sys
import datetime
import os
import chilin

# alias
env = Environment(loader=PackageLoader('chilin', 'template'))
j = os.path.join
e = os.path.exists
error = logging.error
warn = logging.warning


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
                            temp.append(self.outputd + '/' + self.Nameconfigs[session][option].replace('${control_rep}', str(control_rep)))
                        self.Nameconfigs[session][option] = temp

                if self.treat_path[0] != '':
                    if '${treat_rep}' in self.Nameconfigs[session][option]:
                        for treat_rep in range(1, len(self.treat_path) + 1):
                            temp.append(self.Nameconfigs[session][option].replace('${treat_rep}', str(treat_rep)))
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

# parent class of DC part
class PipeController(object):

    def __init__(self):
        """
        read in Options from command line
        Get template and conf information
        """
        self.has_run = True
        self.cmd = ''
        self.shellextract = ''

    def run(self):
        self.shellextract = call(self.cmd, shell = True)
        if self.shellextract !=0:
            self.has_run = False

    def partition(self, step, target = ''):
        '''use step as folder'''
        if self.has_run:
            if not e():
                self.cmd = 'mkdir %s' % (step)
            else:
                self.cmd = 'mv %s %s' %(target, step)
            self.run()
            mv = call(self.cmd, shell = True)
            if mv != 0:
                logging.error('partition error')
                self.has_run = False

    def render(self, template = ''):
        """
        write into the DA.txt template
        """
        if self.has_run:
            print 'write the template'

# children class
class PipeBowtie(PipeController):
    """Bowtie DC and QC step"""
    def __init__(self, chilinconfigs, nameconfigs):
        super(PipeBowtie, self).__init__()
        self.chilinconfigs = chilinconfigs
        self.nameconfigs = nameconfigs
        self.bowtieinfo = {}

    def _format(self, files):
        '''format input and output format'''

        for rep in range(len(files)):
            self.cmd = '{0} view -bt {1} {2} -o {3}'
            self.cmd = self.cmd.format(self.chilinconfigs['samtools']['samtools_main'],
                                       self.chilinconfigs['samtools']['samtools_chrom_len_path']+\
                                           'chromInfo_' + self.chilinconfigs['userinfo']['species'] + '.txt',
                                       files[rep],
                                       self.nameconfigs['bowtieresult']['bam_treat'],
                                       )
            self.run()

    def extract(self):
        '''
        extract information from sam file
        to write data summary and qc measurement
        revision from Xikun and Junsheng
        '''

        for sam_rep in xrange(len(self.nameconfigs['bowtietmp']['treat_sam'])):
            samfile = file(self.nameconfigs['bowtietmp']['treat_sam'][sam_rep],"r")
            reads_dict = {}
            location_dict = {}
            total_reads = 0
            mapped_reads = 0

            for line in samfile:
                sam_line = line.split("\t")
                if sam_line[0] in ['@HD','@SQ','@PG','@RG']: # eliminate the table's head
                    continue
                else:
                    total_reads += 1
                    if sam_line[1] == "4": # "4" means not mapped
                        continue
                    else:
                        mapped_reads += 1 # mapped reads
                        #location = sam_line[1]+sam_line[2]+sam_line[3]
                location = sam_line[1]+":"+sam_line[2]+":"+sam_line[3] ####edited by bo to avoid mistakes20111228
                read_name = sam_line[0]

                if reads_dict.has_key(read_name):
                    reads_dict[read_name] += 1
                elif not reads_dict.has_key(read_name):
                    reads_dict[read_name] = 1
                else:
                    print("! STRANGE in reads_dict")

                if location_dict.has_key(location):
                    location_dict[location] += 1
                elif not location_dict.has_key(location):
                    location_dict[location] = 1
                else:
                    print("! STRANGE in location_dict")

            uniq_read = 0
            uniq_location = 0
            for read_name in reads_dict.keys():####edited by bo change == to >=
                if reads_dict[read_name] >= 1: # the read was only mapped once
                    uniq_read += 1
            for location in location_dict.keys():
                if location_dict[location] >= 1: # the location was only covered once
                    uniq_location += 1
            usable_percentage = float(uniq_read)/float(total_reads)*100
        self.bowtieinfo['mapped'] = mapped_reads
        self.bowtieinfo['unique'] = uniq_read
        self.bowtieinfo['uniq_location'] = uniq_location
        self.bowtieinfo['uniq_ratio'] = str(usable_percentage) + '%'
        self.bowtieinfo['total'] = total_reads

    def process(self):
        treatpath = self.chilinconfigs['userinfo']['treatpath'].split(',')
        for treat_rep in range(len(treatpath)):
            self.cmd  = '{0} -S -m {1} {2} {3} {4} '
            self.cmd = self.cmd.format(self.chilinconfigs['bowtie']['bowtie_main'],
                                        self.chilinconfigs['bowtie']['nbowtie_max_alignment'],
                                        self.chilinconfigs['bowtie']['bowtie_genome_index_path'] + \
                                            self.chilinconfigs['userinfo']['species'],
                                        treatpath[treat_rep],
                                        self.nameconfigs['bowtietmp']['treat_sam'][treat_rep])
            print "bowtie is processing %s" %(treatpath[treat_rep])
            self.run()
            self._format(self.nameconfigs['bowtietmp']['treat_sam'])
            self.extract()
        self.render('test')



class PipeMACS2(PipeController):
    """ MACS step, separately and merge for sorted bam
    for peaks calling
    macs2 callpeak -B -q 0.01 --keep-dup 1 --shiftsize=73 --nomodel  -t /Users/qianqin/Documents/testchilin/testid_treat_rep2.sam -n test.bed"""

    def __init__(self):
        super(PipeMACS2, self).__init__()

    def _format(self):
        return

    def process(self,):
        cmd = '{0} callpeak -B -q 0.01 --keep-dup 1 --shift-size {2} --nomodel ' + \
              '-t {3} {4} -n {5}'


        cmd = cmd.format(self.macs2_main,
                         self.genome_option,
                         self.shiftsize,
                         self.treat_bam,
                         self.control_bam,
                         self.macsname)

        return cmd


class PipeVennCor(PipeController):
    def __init__(self, OptionMethod = "Mean"):
        """Read in the Controller Dictionary
        to decide whether do this step or not"""
        super(PipeVennCor, self).__init__()
        print "if replicates, Do this Step"

    def _format(self, peaksnumber):
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
