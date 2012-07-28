from subprocess import call
from jinja2 import Environment, PackageLoader
from ConfigParser import ConfigParser
import logging
import string
import datetime
import os
import chilin
import re

# alias
env = Environment(loader=PackageLoader('chilin', 'template'))
j = os.path.join
e = os.path.exists
error = logging.error
warn = logging.warning

def repanalysis(datafilelist, filetype = ''):
    '''to reduce code'''
    for rep in datafilelist:
        print 'treat or control'

    print 'analysis treats or controls in batch'


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
        print self.ChiLinconfigs
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
        self.cf = ConfigParser()
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
        self.logger = logging.getLogger()
        self.logfile = logfile
        handler = logging.FileHandler(logfile)
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s :  ')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.NOTSET)

    def record(self, logcontent):
#        dt = datetime.datetime.now()
        self.logger.info(logcontent)
#        with open(self.logfile, 'wb') as f:
#            f.write(dt.strftime('%Y-%m-%d-%H:%M:%S') + logcontent)

# parent class of DC part
class PipeController(object):
    def __init__(self):
        """
        read in Options from command line
        Get template and conf information
        """
        self.has_run = True
        self.cmd = ''
        self.shellextract = 0

    def run(self):
        self.shellextract = call(self.cmd, shell = True)
        if self.shellextract !=0:
            self.has_run = False
        else:
            self.has_run = True

    def partition(self, step, target = '', newname = ''):
        '''use step as folder'''
        if self.has_run:
            if not e(step):
                self.cmd = 'mkdir %s' % (step)
                self.run()
                self.cmd = 'mv %s %s' %(target, step + '/' + newname)
                self.run()
            else:
                self.cmd = 'mv %s %s' %(target, step + '/' + newname)
                self.run()
        print self.cmd

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
        self.treatpath = ''
        self.bowtieinfo = {}

    def _format(self, sam, bam):
        '''format input and output format'''
        if self.has_run:
            self.cmd = '{0} view -bt {1} {2} -o {3}'
            self.cmd = self.cmd.format(self.chilinconfigs['samtools']['samtools_main'],
                                       self.chilinconfigs['samtools']['samtools_chrom_len_path']+\
                                                'chromInfo_' + self.chilinconfigs['userinfo']['species'] + '.txt',
                                       sam,
                                       bam,
                                       )
            self.run()

    def extract(self):
        '''
        extract information from sam file
        to write data summary and qc measurement
        revision from Xikun and Junsheng
        '''
        for sam_rep in range(len(self.nameconfigs['bowtietmp']['treat_sam'])):
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
            self.bowtieinfo['treatrep%s_mapped'%(sam_rep + 1)] = mapped_reads
            self.bowtieinfo['treatrep%s_unique'%(sam_rep + 1)] = uniq_read
            self.bowtieinfo['treatrep%s_uniq_location'%(sam_rep + 1)] = uniq_location
            self.bowtieinfo['treatrep%s_uniq_ratio'%(sam_rep + 1)] = str(usable_percentage) + '%'
            self.bowtieinfo['treatrep%s_total'%(sam_rep + 1)] = total_reads

    def process(self):
        treatpath = self.chilinconfigs['userinfo']['treatpath'].split(',')
        os.system('mkdir bowtie')
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


            self._format(self.nameconfigs['bowtietmp']['treat_sam'][treat_rep], self.nameconfigs['bowtieresult']['bam_treat'][treat_rep])
            self.partition('bowtie', self.nameconfigs['bowtietmp']['treat_sam'][treat_rep])
            self.partition('bowtie', self.nameconfigs['bowtieresult']['bam_treat'][treat_rep])
        os.chdir('bowtie')
        self.extract()
        os.chdir('..')
        self.render()

class PipeMACS2(PipeController):
    """ MACS step, separately and merge for sorted bam
    for peaks calling
    macs2 callpeak -B -q 0.01 --keep-dup 1 --shiftsize=73 --nomodel  -t /Users/qianqin/Documents/testchilin/testid_treat_rep2.sam -n test.bed"""

    def __init__(self, chilinconfigs, nameconfigs):
        super(PipeMACS2, self).__init__()
        self.chilinconfigs = chilinconfigs
        self.nameconfigs = nameconfigs
        self.macsinfo = {}

    def _format(self, bdg, tmp, bw):
        """
        merge bam files
        filter bdg file and remove over-border coordinates
        convert bdg to bw
        command_line = configs["macs.bedgraphtobigwig_main"]+" "+configs["macs.output_treat_bdg_replicates"][i-1]+".tmp "+configs["macs.bedgraphtobigwig_chrom_len_path"]+" "+configs["macs.                    output_treat_bw_replicates"][i-1]
        """
        if self.has_run:
            self.cmd = '{0} -a {1} -b {2} -wa -f 1.00 > {3}'  # bdg filter
            self.cmd = self.cmd.format(self.chilinconfigs['bedtools']['intersectbed_main'],
                                       bdg,
                                       self.chilinconfigs['samtools']['chrom_len_bed_path'] + self.chilinconfigs['userinfo']['species'] + '.bed',
                                       tmp)
            self.run()
        if self.has_run:
            self.cmd = '{0} {1} {2} {3} ' # bedGraphTobigwiggle
            self.cmd = self.cmd.format(self.chilinconfigs['macs']['bedgraphtobigwig_main'],
                                       tmp,
                                       self.chilinconfigs['samtools']['samtools_chrom_len_path'] + self.chilinconfigs['userinfo']['species'] + '.txt', # bedGraphTobw chrom len equal to samtools'
                                       bw,
                                       )
            self.run()

    def _median(self, nums):
        p = sorted(nums)
        l = len(p)
        if l%2 == 0:
            return (p[l/2]+p[l/2-1])/2
        else:
            return p[l/2]

    def extract(self, peak_file):
        '''
        get macsinfo of peak calling
        from peaks.xls
        add fc_10 and ratio'''
        fhd = open( peak_file,"r" )
        d_wo_FDR = []
        d_w_FDR = []
        t = ""
        for i in fhd:
            i = i.strip()
            if i and not i.startswith("#") and not i.startswith("chr\t"):
                fs = i.split("\t")
                fc = fs[7]
                if len(fs) == 9:
                    # has FDR
                    d_w_FDR.append((float(fc),float(fs[8])))
                else:
                    d_wo_FDR.append((float(fc)))
                    if d_w_FDR:
                        d = sorted(d_w_FDR,cmp=lambda x,y:cmp(x[0],y[0]))
                        fdr = [x[1] for x in d if x[0] >= 20] # debug site ####JunSheng debug Sheng'en edit in 20110929, for correct fdr>20%number in summary file
                        fdr_all = [x[1] for x in d]##### Sheng'en edit in 20111128, collect the condition no fc >= 20
                        t += "peaks_fc_ge_20 = %d\n" % len(fdr)
                        if len(fdr)>0:
                            t += "max_FDR_of_fc_ge_20 = %.2f\n" % max(fdr)
                            t += "median_FDR_of_fc_ge_20 = %.2f\n" % median(fdr)
                        else:
                            t += "No peaks FDR fc>=20\nmax_FDR_of_fc = %.2f\n"% max(fdr_all)
                            t += "median_FDR_of_fc = %.2f\n" % median(fdr_all)
                    else:
                        d = sorted(d_wo_FDR)
                        d = [x for x in d if x>=20]
                        t += "peaks_fc_ge_20 = %d\n" % len(d)####20120514 modified by qinqian, from 'le_20' to 'ge_20'
                        return t

    def process(self, shiftsize = ''):
        # macs2, original generated results
        '''testid_macs_control_lambda.bdg # temp
        testid_macs_peaks.bed # result
        testid_macs_peaks.encodePeak # temp
        testid_macs_peaks.xls  # result
        testid_macs_summits.bed # result
        testid_macs_treat_pileup.bdg # -> bw result'''
        for rep in range(len(self.nameconfigs['bowtieresult']['bam_treat'])):
            self.cmd = '{0} callpeak -B -q 0.01 --keep-dup 1 --shiftsize {1} --nomodel ' + \
                  '-t {2} -n {3}'
            self.cmd = self.cmd.format(self.chilinconfigs['macs']['macs_main'],
                                       shiftsize,
                                       'bowtie/' + self.nameconfigs['bowtieresult']['bam_treat'][rep],
                                        # self.control_bam,
                                       self.nameconfigs['macstmp']['macs_initrep_name'][rep])
            self.run()
            if self.has_run:
                self.partition('macs2', self.nameconfigs['macstmp']['macs_initrep_name'][rep] + '_peaks.bed', \
                                   self.nameconfigs['macsresult']['treatrep_peaks'][rep])
                self.partition('macs2', self.nameconfigs['macstmp']['macs_initrep_name'][rep] + '_peaks.xls', \
                                    self.nameconfigs['macsresult']['peaksrep_xls'][rep])
                self.partition('macs2', self.nameconfigs['macstmp']['macs_initrep_name'][rep] + '_summits.bed', \
                                    self.nameconfigs['macsresult']['rep_summits'][rep])
                self.partition('macs2', self.nameconfigs['macstmp']['macs_initrep_name'][rep] + '_treat_pileup.bdg', \
                                    self.nameconfigs['macstmp']['treatrep_bdg'][rep])
            if self.has_run:
                self._format('macs2/' + self.nameconfigs['macstmp']['treatrep_bdg'][rep], \
                                'macs2/' + self.nameconfigs['macstmp']['treatrep_tmp_bdg'][rep], 'macs2/' + self.nameconfigs['macsresult']['treatrep_bw'][rep])

        if len(self.nameconfigs['bowtieresult']['bam_treat']) > 1:
            self.cmd = '{0} merge {1}  {2}'
            self.cmd = self.cmd.format(self.chilinconfigs['samtools']['samtools_main'],
                                       self.nameconfigs['bowtieresult']['bammerge'],
                                       '  '.join(map(lambda x: 'bowtie/' + x, self.nameconfigs['bowtieresult']['bam_treat']))
                                       )
            self.run()
#            self.cmd = self.cmd.format(self.chilinconfigs['samtools']['samtools_main'],
#                                        self.nameconfigs['bowtieresult']['bam_treat'][rep],
#                                        self.chilinconfigs['samtools']['chrom_len_bed_path'] + self.chilinconfigs['userinfo']['species'] + '.bed',
#                                        self.nameconfigs['macstmp']['treatrep_bdg'][rep])
            if self.has_run:
                self.cmd = '{0} callpeak -B -q 0.01 --keep-dup 1 --shiftsize {1} --nomodel ' + \
                       '-t {2} -n {3}'
                self.cmd = self.cmd.format(self.chilinconfigs['macs']['macs_main'],
                                           shiftsize,
                                           self.nameconfigs['bowtieresult']['bammerge'],
                                           self.nameconfigs['macstmp']['macs_init_mergename'])

                self.run()
                if self.has_run:
                    self.partition('macs2', self.nameconfigs['macstmp']['macs_init_mergename'] + '_peaks.bed',
                                   self.nameconfigs['macsresult']['treat_peaks'])
                    self.partition('macs2', self.nameconfigs['macstmp']['macs_init_mergename'] + '_peaks.xls', 
                                   self.nameconfigs['macsresult']['peaks_xls'])
                    self.partition('macs2', self.nameconfigs['macstmp']['macs_init_mergename'] + '_summits.bed', 
                                   self.nameconfigs['macsresult']['summits'])
                    self.partition('macs2', self.nameconfigs['macstmp']['macs_init_mergename'] + '_treat_pileup.bdg', 
                                   self.nameconfigs['macstmp']['treat_bdg'])
                    print 1
                if self.has_run:
                    self._format('macs2/' + self.nameconfigs['macstmp']['treat_bdg'], 
                                 'macs2/' + self.nameconfigs['macstmp']['treat_bdgtmp'], 
                                 'macs2/' + self.nameconfigs['macsresult']['treat_bw'])
               # for control
            if self.chilinconfigs['userinfo']['controlpath'] != '':
                self.cmd = '{0} callpeak -B -q 0.01 --keep-dup 1 --shiftsize {1} --nomodel ' + \
                    '-t {2} -c {3} -n {4}'

class PipeVennCor(PipeController):
    def __init__(self, chilinconfigs, nameconfigs, peaksnumber = '', OptionMethod = "Mean", replicates = False):
        """Read in the Controller Dictionary
        to decide whether do this step or not"""
        super(PipeVennCor, self).__init__()
        self.peaksnumber = ''
        self.ratio = {}
        print "if replicates, Do this Step"

    def _format(self, a_type):
        self.cmd = '{0} -wa -u -a {1} -b {2} > {3}' # intersected
        if a_type == 'dhs':
            self.cmd = self.cmd.format(self.chilinconfigs['bedtools']['intersectbed_main'],
                                       self.nameconfigs['macsresult']['treat_peaks'],
                                       self.chilinconfigs['venn']['dhs_bed_path'],
                                       self.nameconfigs['bedtoolstmp']['dhs_bed'])
        if a_type == 'velcro':
            self.cmd = self.cmd.format(self.chilinconfigs['bedtools']['intersectbed_main'],
                                       self.nameconfigs['macsresult']['treat_peaks'],
                                       self.chilinconfigs['venn']['velcro_path'],
                                       self.nameconfigs['bedtoolstmp']['velcro_bed'])
        self.run()

    def extract(self, a_type):
        """
        extract dhs overlap and velcro overlap information
        """
        lena_type = len(open(a_type, 'r').readlines())
        lenall = len(open(self.nameconfigs['macsresult']['treat_bed']))
        self.ratio[a_type] = float(lena_type) / lenall
        a_type.close()

    def process(self, rep):
        
        self.extract('dhs')
        self.extract('velcro')

        if len(self.chilinconfigs['trepn']) < 2:
            self.cmd = 'touch %s' % self.nameconfigs['represult']['ven_png']
            self.run()
        else:
            # venn diagram
            self.cmd = '{0} -t Overlap_of_Replicates {1} {2}'
            self.cmd.format(self.chilinconfigs['venn']['venn_diagram_main'],
                            ' '.join(self.nameconfigs['macsreult']['treatrep_peaks']),
                            ' '.join(map(lambda x: "-l replicate_" + str(x), \
                                             xrange(1, self.chilinconfigs['trepn'] + 1)))
                            )
            self.run()
            # correlation plot
            if self.has_run:
                self.cmd = '{0} {1} {2} {3} {4} {5}'
                self.cmd.format()
                self.render()


class PipeCEAS(PipeController):
    def __init__(self, chilinconfigs, nameconfigs, peaksnumber):
        """Get CEAS dependency info from
        Check Class"""
        super(PipeCEAS, self).__init__()
        self.chilinconfigs = chilinconfigs
        self.naemconfigs = nameconfigs
        self.peaks = peaksnumber

    def _format(self):
        # get peaksnumber ge >= 10 for ceas
        xls = self.nameconfigs['macsresult']['peaks_xls']
        bed = open(self.nameconfigs['ceasresult']['ceasge10_bed'])
        peaksnum = 0
        for line in open(xls, 'rU'):
            if re.search('chr\w', line[0]):
                line = line.strip()
                line = line.split('\t')
                peaksnum += 1
                if peaksnum <= self.peaks and float(line[7]) >= 10:
                    bed_line = line[0] + '\t' + str(int(line[1]) - 1) + '\t' + line[2] + '\t' +\
                        "macs_peaks_" + str(peaksnum) + '\t' + line[8] + '\n'
                    bed.write(bed_line)
        bed.close()
        xls.close()





    def extract(self):
        '''
        get R code for QC measurement'''
        print "Run the Dependency Program for return FindPath string"

    def process(self):
        print "Call private _Run"
        print "Extract shell output"
        print "Write into the template"

class PipeConserv(PipeController):
    def __init__(self, chilinconfigs, nameconfigs):
        super(PipeConserv, self).__init__()

    def _format(self):
        print "Set input Peaks number for conservation Plot"

    def extract(self):
        print "Run the External conservation plot"

    def process(self):
        print "Call private _Run"
        print "Extract shell output" # generate temporary file
        print "Call FindPath DC to return path for QC"
        print "Call FindPath to return path for DC"
        print "Write into the template"
        return # information for passing to TemplateParser

class PipeMotif(PipeController):
    def __init__(self, chilinconfigs, nameconfigs, peaksnumber):
        super(PipeMotif, self).__init__()
        self.peaksnumber =peaksnumber

    def _format(self):
        # generate top peaks from summits BED file
        summitsfile = open(self.nameconfigs['macsresult']['summits'])
        peaks = []
        for i in summitsfile:
            peaks.append( (i,float(i.split()[-1])) )
            top_n = self.peaksnumber
            top_n_summits = map(lambda x:x[0],sorted(peaks, key=lambda x:x[1], reverse=True)[:top_n])

            top_n_summits_file = "top"+str(top_n)+"_summits.bed"
            top_n_summits_fhd = open(top_n_summits_file,"w")
        for i in top_n_summits:
            top_n_summits_fhd.write(i)
        top_n_summits_fhd.close()

        # remove chrM from top n summits
        f=open(top_n_summits_file,"rU")
        temp=open("2.bed","w")
        for i in f:
            i=i.split()
            if i[0]=="chrM":
                continue
            else:
                temp.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\t"+i[3]+"\t"+i[4]+"\n")
                temp.close()

    def extract(self, zip):
        '''extract information for qc part'''
        print 'get zip information'

    def process(self):

        self._format()
        self.cmd = ''
        self.extract()
        self.render()

class PipeGO(PipeController):
    def __init__(self, chilinconfigs, nameconfigs):
        pass
    def _format(self):
        '''Use RegPotential'''

        pass
class API(PipeController):
    '''internet interface'''
    def __init__(self):
        import cgi
        pass
