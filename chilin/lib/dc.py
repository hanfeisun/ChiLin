import os
import re
import logging
import sys
from subprocess import call
from ConfigParser import SafeConfigParser
from glob import glob
from pkg_resources import resource_filename
from jinja2 import Environment, PackageLoader

exists = os.path.exists
error = logging.error
warn = logging.warning

class LogWriter:
    def __init__(self, logfile = 'log'):
        """
        Universal log format
        time + incidence
        """
        self.logger = logging.getLogger()
        self.logfile = logfile
        handler = logging.FileHandler(logfile)
        formatter = logging.Formatter('%(asctime)s %(levelname)s : %(message)s ;  ')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.NOTSET)

    def record(self, logcontent):
        print logcontent
        self.logger.info(logcontent)
        return logcontent

def gen_conf( species ):
    env = Environment(loader = PackageLoader('chilin', 'template'))
    temp = env.get_template('ChiLinjinja.conf')
    if species == 'hg19':
        conf = temp.render(species = species,
                           filterdup_species = 'hs')
    elif species == 'mm9':
        conf = temp.render(species = species,
                           filterdup_species = 'mm')
    print conf

class PipePreparation:
    def __init__(self, ChiLinconfPath,
                 NameConfPath = resource_filename("chilin", os.path.join("db", "NameRule.conf"))):
        """
        Parse the Name Rule and 
        Customer filled Chilinconf
        """
        self.ChiLinconfPath = ChiLinconfPath
        self.NameConfPath = NameConfPath
        self._rule = {}
        self._conf = {}
        self.checked = False
        self._read_conf(self.ChiLinconfPath)
        parseinput = lambda x: x.split(',')
        self._conf['userinfo']['treatpath'] = parseinput(self._conf['userinfo']['treatpath'])
        self._conf['userinfo']['controlpath'] = parseinput(self._conf['userinfo']['controlpath'])
        self._conf['userinfo']['treatnumber'] = len(self._conf['userinfo']['treatpath'])
        self._conf['userinfo']['controlnumber'] = len(self._conf['userinfo']['controlpath'])
        os.chdir(self._conf['userinfo']['outputdirectory'])

        self._read_rule(self.NameConfPath)
        self.log = LogWriter(self._rule['root']['log']).record

    def get_config(self):
        if self.checked == False:
            self.check_conf()
        return self._conf
        
    def get_rule(self):
        if self.checked == False:
            self.check_conf()
        return self._rule
        
    def get_tex(self):
        return self._rule["qcresult"]["qctex"]

    def get_summary(self):
        return self._rule['root']['data_summary']
        
    def func_log(self):
        return self.log
        
    def _read_rule(self, rule_path):
        """
        Read in the conf of NameRule.Conf with replacement of %(DatasetID), %(treat_rep) and %(control_rep)
        """

        cf = SafeConfigParser()
        cf.read(rule_path)
        cf.set('DEFAULT',
               'DatasetID', self._conf['userinfo']['datasetid'])
        uni = lambda str: str.strip().lower()
        
        for sec in cf.sections():
            get_raw = lambda opt: cf.get(sec, opt, 1)
            get_expand = lambda opt: cf.get(sec, opt, 0)
            def get_fmt(opt):
                get_rep_expand = lambda raw_str, rep_cnt: map(lambda x:cf.get(sec, opt, 0, {raw_str: str(x+1)}),
                                                              range(rep_cnt))
                if '(treat_rep)' in get_raw(opt):
                    return get_rep_expand("treat_rep",
                                          self._conf['userinfo']['treatnumber'])
                elif '(control_rep)' in get_raw(opt):
                    return get_rep_expand("control_rep",
                                          self._conf['userinfo']['controlnumber'])
                else:
                    return get_expand(opt)
                    
            self._rule[uni(sec)] = dict((uni(opt), get_fmt(opt)) for opt in cf.options(sec))
        
    def _read_conf(self, conf_path):
        """
        Read in ChiLin.conf without replacement
        """
        cf = SafeConfigParser() 
        cf.read(conf_path)
        uni = lambda str: str.strip().lower()
        for sec in cf.sections():
            get_raw = lambda opt: cf.get(sec, opt, 1)
            self._conf[uni(sec)] = dict((uni(opt), get_raw(opt)) for opt in cf.options(sec))


    def check_conf(self):
        """
        Check the Meta configuration
        if up to our definition
        """

        # Convert treat, control in to list
        # get replicates number
        head = lambda a_list: a_list[0] != ""

        def fatal(prediction, error_info,
                  side_effect=lambda : True):
            def check():
                if prediction:
                    error(error_info)
                    sys.exit(1)
                else:
                    side_effect()
            return check
        def danger(prediction, warn_info):
            def check():
                if prediction:
                    warn(warn_info)
            return check
        
        check_list = lambda check_funcs: (f() for f in check_funcs)
        
        c = []
        c.append(fatal(not head(self._conf['userinfo']['treatpath']),
                       'forget to fill the treat file path'))
        c.append(danger(not head(self._conf['userinfo']['controlpath']),
                        'No control file path'))
        c.append(fatal(not os.path.isdir(self._conf['userinfo']['outputdirectory']),
                       'check output directory name'))
        c.append(fatal(self._conf['userinfo']['species'] not in ['hg19', 'mm9'],
                       'forget to fill the species or species input not supported'))
        c.append(fatal(False in map(os.path.isfile, self._conf['userinfo']['treatpath']),
                       'check your input treat file, some error'))
        c.append(fatal(False in map(os.path.isfile, self._conf['userinfo']['controlpath']),
                       'check your input control file, some error'))
        c.append(fatal(not exists(self._conf['qc']['fastqc_main']),
                       'fastqc not exists'))
        c.append(fatal(not exists(self._conf['bowtie']['bowtie_main']),
                       "bowtie program dependency has problem"))
        c.append(fatal(not exists(self._conf['samtools']['samtools_main']),
                       "samtools dependency has problem"))
        c.append(fatal(not exists(self._conf['macs']['macs_main']),
                       "macs2 program dependency has probelm"))
        c.append(fatal(not exists(self._conf['bedtools']['intersectbed_main']),
                       "bedtools program dependency has problem"))
        c.append(fatal(not exists(self._conf['conservation']['conserv_plot_main']),
                       "conservation_plot dependency has problem"))
        c.append(fatal(not exists(self._conf['correlation']['wig_correlation_main']),
                       "correlation plot dependency has problem"))
        c.append(fatal(not exists(self._conf['venn']['venn_diagram_main']),
                       "venn plot dependency has problem"),)

        check_list(c)
        self.log('ChiLin config parse success, and dependency check has passed!')        
        self.checked = True


class PipeController(object):
    def __init__(self, conf, rule, log, **args):
        """
        read in Options from command line
        Get template and conf information
        """
        self.cmd = ''
        self.shellextract = 0
        self.env = Environment(loader=PackageLoader('chilin', 'template'))
        self.conf = conf
        self.rule = rule
        self.log = log
        self.debug = args.get("debug", False)
        self.threads = args.get("threads", 1)
        self.datasummary = args.get("datasummary", 1)

    def run_cmd(self, cmd, exit_ = True, error_handler = lambda :True):
        """
        univeral call shell and judge
        """
        self.log("Run command:\t"+cmd)
        if call(cmd, shell = True):
            # if encounters error
            result = error_handler()
            if exit_:
                print "`"+cmd+"`"+" failed, exit"
                sys.exit(0)
            else:
                return result
        else:
            return True

    def if_runcmd(self, condition, cmd, else_handler = lambda :True, size_check=True):
        """
        run a command conditionally
        """
        if type(condition) == str:
            self.log("checking file "+condition)
            if not exists(condition):
                condition = True
            else:
                if os.path.isfile(condition) and os.path.getsize(condition) <=0 and size_check:
                    condition = True
                else:
                    condition = False
                    
        if condition:
            return self.run_cmd(cmd)
        else:
            else_handler()
            return self.log(cmd+" is skipped")


        
    def _render(self):
        """
        write into the DA.txt template
        """
        DA = self.env.get_template('DA.txt')
        ds_rendered = DA.render(self.rendercontent)
        with open(self.datasummary,"a") as df:
            df.write(ds_rendered)
class PipeGroom(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, **args):
        """
        pipeline bowtie part"""
        super(PipeGroom, self).__init__(conf, rule, log, **args)
        self.stepcontrol = stepcontrol
    def run(self):
        groom_path = lambda x:x.replace(".bam", ".fastq")
        need_groom = lambda x:".bam" in x
        for treat_rep in range(self.conf['userinfo']['treatnumber']):
            if need_groom(self.conf['userinfo']['treatpath'][treat_rep]):
                cmd  = '{0} -q {1} > {2}'
                cmd = cmd.format(self.conf['bowtie']['BAMTOFASTQ'],
                                 self.conf['userinfo']['treatpath'][treat_rep],
                                 groom_path(self.conf['userinfo']['treatpath'][treat_rep]))
                self.log("bam2fastq is processing %s" % (self.conf['userinfo']['treatpath'][treat_rep]))
                self.run_cmd(cmd)
                self.conf['userinfo']['treatpath'][treat_rep] = groom_path(self.conf['userinfo']['treatpath'][treat_rep])
        for control_rep in range(self.conf['userinfo']['controlnumber']):
            if need_groom(self.conf['userinfo']['controlpath'][control_rep]):
                cmd  = '{0} -q {1} > {2}'
                cmd = cmd.format(self.conf['bowtie']['BAMTOFASTQ'],
                                 self.conf['userinfo']['controlpath'][control_rep],
                                 groom_path(self.conf['userinfo']['controlpath'][control_rep]))
                self.log("bam2fastq is processing %s" % (self.conf['userinfo']['controlpath'][control_rep]))
                self.run_cmd(cmd)
                self.conf['userinfo']['controlpath'][control_rep] = groom_path(self.conf['userinfo']['controlpath'][control_rep])
        


class PipeBowtie(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, **args):
        """
        pipeline bowtie part"""
        super(PipeBowtie, self).__init__(conf, rule, log, **args)
        self.rendercontent = {}
        self.stepcontrol = stepcontrol
        
    def _sam2bam(self, sam, bam):
        """
        using samtools
        convert bowtie sam result to bam
        example: samtools view -bt chrom_len sam bam
        """
        cmd = '{0} view -bt {1} {2} -o {3}'
        cmd = cmd.format(self.conf['samtools']['samtools_main'],                          
                         self.conf['samtools']['samtools_chrom_len_path'], # chrom_len depends on species
                         sam, bam)
        if self.debug:
            self.if_runcmd(bam, cmd)
        else:
            self.run_cmd(cmd)
                           
    

    def _extract(self, cnt, files, control = False):
        """
        inputfile : sam files(treat or control)
        cnt : replicates number
        extract bowtie qc information
        sams = [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...] **args
        sams = self.rendercontent
        sams = {'sams': [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...]} *arg
        """
        print "working on extracting"
        for sam_rep in range(cnt):
            reads_dict = {}
            location_dict = {}
            total_reads = 0
            mapped_reads = 0
            with open(files[sam_rep]) as samfile:
                for line in samfile:
                    if line.startswith("@"): # eliminate the table's head
                        continue
                    li = line.split("\t")
                    
                    total_reads += 1
                    if total_reads % 1000000 == 0: print total_reads
                    if len(li) < 3:
                        print line + "looks strange, skipped"
                        continue
                    if li[1] == "4": # "4" means not mapped
                        continue

                    
                    mapped_reads += 1
                    location = ":".join(li[1:4])
                    reads_dict[li[0]] = reads_dict.get(li[0], 0) + 1
                    location_dict[location] = location_dict.get(location, 0) + 1
                    
            uniq_read = len([i for i in reads_dict if reads_dict[i] == 1])
            uniq_location = len(location_dict)
            usable_percentage = float(uniq_read)/float(total_reads)*100

            info = { 'name':self.rule['bowtietmp']['treat_sam'][sam_rep] if not control else self.rule['bowtietmp']['control_sam'][sam_rep],
                     'total': total_reads,
                     'mapped': mapped_reads,
                     'unireads': uniq_read,
                     'uniloc': uniq_location,
                     'percentage' : str(usable_percentage) + '%'}
            
            self.rendercontent['sams'].append(info)

    def run(self):
        self.rendercontent['sams'] = []
        for treat_rep in range(self.conf['userinfo']['treatnumber']):
            cmd  = '{0} -p {5} -S -m {1} {2} {3} {4} '
            cmd = cmd.format(self.conf['bowtie']['bowtie_main'],
                             self.conf['bowtie']['nbowtie_max_alignment'],
                             self.conf['bowtie']['bowtie_genome_index_path'],
                             self.conf['userinfo']['treatpath'][treat_rep],
                             self.rule['bowtietmp']['treat_sam'][treat_rep],
                             self.threads)
            self.log("bowtie is processing %s" % (self.conf['userinfo']['treatpath'][treat_rep]))
            if self.debug:
                self.if_runcmd(self.rule['bowtietmp']['treat_sam'][treat_rep], cmd)
            else:
                self.run_cmd(cmd)
            self._sam2bam(self.rule['bowtietmp']['treat_sam'][treat_rep], self.rule['bowtieresult']['bam_treat'][treat_rep])


        for control_rep in range(self.conf['userinfo']['controlnumber']):
            cmd  = '{0} -p {5} -S -m {1} {2} {3} {4} '
            cmd = cmd.format(self.conf['bowtie']['bowtie_main'],
                             self.conf['bowtie']['nbowtie_max_alignment'],
                             self.conf['bowtie']['bowtie_genome_index_path'],
                             self.conf['userinfo']['controlpath'][control_rep],
                             self.rule['bowtietmp']['control_sam'][control_rep],
                             self.threads)
            self.log("bowtie is processing %s" % (self.conf['userinfo']['controlpath'][control_rep]))
            if self.debug:
                self.if_runcmd(self.rule['bowtietmp']['control_sam'][control_rep], cmd)
            else:
                self.run_cmd(cmd)
            self._sam2bam(self.rule['bowtietmp']['control_sam'][control_rep], self.rule['bowtieresult']['bam_control'][control_rep])
        if exists(self.datasummary) and self.debug:
            print "skip rendering"
            pass
        else:
            self._extract(self.conf['userinfo']['treatnumber'], self.rule['bowtietmp']['treat_sam'], control = False)
            self._extract(self.conf['userinfo']['controlnumber'], self.rule['bowtietmp']['control_sam'], control = True)
            self._render()
        self.log("bowtie run successfully")

class PipeMACS2(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, shiftsize, **args):
        """
        MACS step, separately and merge for sorted bam
        shell example:
        macs2 callpeak -B -q 0.01 --keep-dup 1 \
                --shiftsize=73 --nomodel  -t /Users/qianqin/Documents/testchilin/testid_treat_rep2.sam  -c control.bam -n test.bed
        """
        super(PipeMACS2, self).__init__(conf, rule, log, **args)
        self.macsinfo = {}

        self.stepcontrol = stepcontrol
        self.shiftsize = shiftsize
        self.rendercontent = {}

    def _format(self, bdg, tmp, bw):
        """
        filter bdg file and remove over-border coordinates
        convert bdg to bw
        shell example
        /usr/local/bin/intersectBed -a 6523_rep1_treat.bdg \ 
                -b /opt/bin/chr_limit/chr_limit_mm9.bed -wa -f 1.00 > 6523_rep1_treat.bdg.tm
        /opt/bin/UCSCTools/bedGraphToBigWig 6523_control.bdg.tmp \
                /mnt/Storage/data/Samtool/chromInfo_mm9.txt 6523_control.bw 
        bedGraphTobw chrom len equal to samtools'
        """
        cmd = '{0} -a {1} -b {2} -wa -f 1.00 > {3}'  # bdg filter
        cmd = cmd.format(self.conf['bedtools']['intersectbed_main'],
                         bdg,
                         self.conf['samtools']['chrom_len_bed_path'],
                         tmp)
        if self.debug:
            self.if_runcmd(tmp, cmd)
        else:
            self.run_cmd(cmd)
        cmd = '{0} {1} {2} {3} ' # bedGraphTobigwiggle
        cmd = cmd.format(self.conf['macs']['bedgraphtobigwig_main'],
                         tmp,
                         self.conf['samtools']['samtools_chrom_len_path'],
                         bw )
        if self.debug:
            self.if_runcmd(bw, cmd)
        else:
            self.run_cmd(cmd)

    def extract(self):
        """
        get macsinfo of peak calling
        from peaks.xls
        add fc_10 and ratio
        macs2info = {'ratios': {'totalpeak':a, 'total1': 5...}..} *args
        """
        fhd = open(self.rule['macsresult']['peaks_xls'],"r")
        total = 0
        fc20n = 0
        fc10n = 0
        for i in fhd:
            i = i.strip()
            if i and not i.startswith("#") and not i.startswith("chr\t"):
                total += 1
                fs = i.split("\t")
                fc = fs[7]
                if fc >= 20:
                    fc20n += 1
                if fc >= 10:
                    fc10n += 1
        self.macsinfo['totalpeak'] = total
        self.macsinfo['peaksge20'] = fc20n
        self.macsinfo['peaksge10'] = fc10n
        if total != 0:
            self.macsinfo['peaksge10ratio'] = float(fc10n)/total
            self.macsinfo['peaksge20ratio'] = float(fc20n)/total
        else:
            self.macsinfo['peaksge10ratio'] = 0.0001
            self.macsinfo['peaksge20ratio'] = 0.0001
        self.macsinfo['distance'] = float(self.shiftsize) * 2
        self.rendercontent['ratios'] = self.macsinfo

    def run(self):
        """
        testid_macs_control_lambda.bdg # temp
        testid_macs_peaks.encodePeak   # temp
        testid_macs_peaks.bed        # result
        testid_macs_peaks.xls        # result
        testid_macs_summits.bed      # result
        testid_macs_treat_pileup.bdg # -> bw result
        shell example
            macs2 callpeak -g hs -B -q 0.01 --keep-dup 1 --nomodel --shiftsize=73 -t GSM486702_BI.CD34_Primary_Cells.Input.CD34_39661.bed -c GSM486702_BI.CD34_Primary_Cells.Input.CD34_39661.bed  -n test1
        """
        if self.stepcontrol < 2:
            sys.exit(1)
        # genome option
        if 'hg' in self.conf['userinfo']['species']:
            genome_option = ' -g hs '
        elif 'mm' in self.conf['userinfo']['species']:
            genome_option = ' -g mm '
        else:
            genome_option = ' '
        if self.conf['userinfo']['treatnumber'] > 1:
            cmd = '{0} merge -f {1}  {2}'
            cmd = cmd.format(self.conf['samtools']['samtools_main'],
                                       self.rule['bowtieresult']['bamtreatmerge'],
                                       '  '.join(self.rule['bowtieresult']['bam_treat']))
            if self.debug:
                self.if_runcmd(self.rule['bowtieresult']['bamtreatmerge'], cmd)
            else:
                self.run_cmd(cmd)

        if self.conf['userinfo']['controlnumber'] > 1:
            cmd = '{0} merge -f {1}  {2}'
            cmd = cmd.format(self.conf['samtools']['samtools_main'],
                                       self.rule['bowtieresult']['bamcontrolmerge'],
                                       '  '.join(self.rule['bowtieresult']['bam_control']))
            if self.debug:
                self.if_runcmd(self.rule['bowtieresult']['bamcontrolmerge'], cmd)
            else:
                self.run_cmd(cmd)

        # set control option, merge bam control for universal control
        if self.conf['userinfo']['controlnumber'] == 1:
            control_option = ' -c ' + self.rule['bowtieresult']['bam_control'][0]
        elif self.conf['userinfo']['controlnumber'] > 1:
            control_option = ' -c ' + self.rule['bowtieresult']['bamcontrolmerge']
        else:
            control_option = ' '

        # each bam file peak calling
        for treat_rep in range(self.conf['userinfo']['treatnumber']):
            cmd = '{0} callpeak {1}  -B -q 0.01 --keep-dup 1 --shiftsize {2} --nomodel ' + \
                  '-t {3} {4} -n {5}'
            cmd = cmd.format(self.conf['macs']['macs_main'],
                                       genome_option,
                                       self.shiftsize,
                                       self.rule['bowtieresult']['bam_treat'][treat_rep],
                                       control_option,
                                       self.rule['macstmp']['macs_initrep_name'][treat_rep])
            # convert macs default name to NameRule
            if self.debug:
                self.if_runcmd(not any([exists(self.rule['macstmp']['macs_initrep_name'][treat_rep]+ '_treat_pileup.bdg'),
                                        exists(self.rule['macstmp']['treatrep_bdg'][treat_rep])]), cmd)
            else:
                self.run_cmd(cmd)
            cmd = 'mv %s %s' % (self.rule['macstmp']['macs_initrep_name'][treat_rep]+ '_treat_pileup.bdg',
                                self.rule['macstmp']['treatrep_bdg'][treat_rep])
            if self.debug:
                self.if_runcmd(self.rule['macstmp']['treatrep_bdg'][treat_rep], cmd)
            else:
                self.run_cmd(cmd)
                
            cmd = 'mv %s %s' % (self.rule['macstmp']['macs_initrep_name'][treat_rep] + '_control_lambda.bdg',
                                self.rule['macstmp']['controlrep_bdg'][treat_rep])
            if self.debug:
                self.if_runcmd(self.rule['macstmp']['controlrep_bdg'][treat_rep], cmd)
            else:
                self.run_cmd(cmd)
            self._format(self.rule['macstmp']['treatrep_bdg'][treat_rep], 
                         self.rule['macstmp']['treatrep_tmp_bdg'][treat_rep], 
                         self.rule['macsresult']['treatrep_bw'][treat_rep])
            self._format(self.rule['macstmp']['controlrep_bdg'][treat_rep],
                         self.rule['macstmp']['controlrep_tmp_bdg'][treat_rep], 
                         self.rule['macsresult']['controlrep_bw'][treat_rep])

        # merge treat bam files for peaks calling
        if len(self.rule['bowtieresult']['bam_treat']) > 1:
            cmd = '{0} callpeak {1} -B -q 0.01 --keep-dup 1 --shiftsize {2} --nomodel -t {3} {4} -n {5}'
            cmd = cmd.format(self.conf['macs']['macs_main'],
                                       genome_option,
                                       self.shiftsize,
                                       self.rule['bowtieresult']['bamtreatmerge'],
                                       control_option,
                                       self.rule['macstmp']['macs_init_mergename'])
            if self.debug:
                self.if_runcmd(not any([exists(self.rule['macstmp']['macs_init_mergename'][treat_rep]+ '_treat_pileup.bdg'),
                                        exists(self.rule['macstmp']['treat_bdg'][treat_rep])]), cmd)
            else:
                self.run_cmd(cmd)
            # convert macs default name to NameRule
            cmd = 'mv %s %s' % (self.rule['macstmp']['macs_init_mergename'] + '_treat_pileup.bdg',
                                self.rule['macstmp']['treat_bdg'])
            if self.debug:
                self.if_runcmd(self.rule['macstmp']['treat_bdg'], cmd)
            else:
                self.run_cmd(cmd)
                
            cmd = 'mv %s %s' % (self.rule['macstmp']['macs_init_mergename'] + '_control_lambda.bdg',
                                self.rule['macstmp']['control_bdg'])
            if self.debug:
                self.if_runcmd( self.rule['macstmp']['control_bdg'], cmd)
            else:
                self.run_cmd(cmd)
                
            self._format(self.rule['macstmp']['treat_bdg'], 
                         self.rule['macstmp']['treat_bdgtmp'], 
                         self.rule['macsresult']['treat_bw'])
            self._format(self.rule['macstmp']['control_bdg'], 
                         self.rule['macstmp']['control_tmp_bdg'], 
                         self.rule['macsresult']['control_bw'])
        elif len(self.rule['bowtieresult']['bam_treat']) == 1:
                cmd = 'cp %s %s' % (self.rule['macsresult']['peaksrep_xls'][0],
                                    self.rule['macsresult']['peaks_xls'])
                self.run_cmd(cmd)
                cmd = 'cp %s %s' % (self.rule['macsresult']['treatrep_peaks'][0],
                                    self.rule['macsresult']['treat_peaks'])
                self.run_cmd(cmd)
                cmd = 'cp %s %s' % (self.rule['macsresult']['rep_summits'][0],
                                    self.rule['macsresult']['summits'])
                self.run_cmd(cmd)
                cmd = 'cp %s %s' % (self.rule['macsresult']['rep_summits'][0],
                                    self.rule['macsresult']['summits'])

                cmd = 'cp %s %s' % (self.rule['macsresult']['treatrep_bw'][0],
                                    self.rule['macsresult']['treat_bw'])
                if self.debug:
                    self.if_runcmd(self.rule['macsresult']['treat_bw'], cmd)
                else:
                    self.run_cmd(cmd)
                    
                cmd = 'cp %s %s' % (self.rule['macsresult']['controlrep_bw'][0],
                                    self.rule['macsresult']['control_bw'])
                if self.debug:
                    self.if_runcmd(self.rule['macsresult']['control_bw'], cmd)
                else:
                    self.run_cmd(cmd)
                
        self.extract()

class PipeVennCor(PipeController):
    def __init__(self, conf, rule, log,
                 stepcontrol, ratios, peaksnumber = '', OptionMethod = "Mean", **args):
        """
        replicates qc measurement 
        using venn diagram and correlation 
        """
        super(PipeVennCor, self).__init__(conf, rule, log, **args)
        self.peaksnumber = peaksnumber
        self.OptionMethod = OptionMethod
        self.ratio = {}
        self.stepcontrol = stepcontrol
        self.rendercontent = ratios

    def _format(self, a_type):
        """
        input : peaks.bed

        get intersect regions between merge peaks bed and velcro, DHS sites
        """
        cmd = '{0} -wa -u -a {1} -b {2} > {3}' # intersected
        if a_type == 'dhs':
            cmd = cmd.format(self.conf['bedtools']['intersectbed_main'],
                                       self.rule['macsresult']['treat_peaks'],
                                       self.conf['venn']['dhs_bed_path'],
                                       self.rule['bedtoolstmp']['dhs_bed'])
        if a_type == 'velcro':
            cmd = cmd.format(self.conf['bedtools']['intersectbed_main'],
                                       self.rule['macsresult']['treat_peaks'],
                                       self.conf['venn']['velcro_path'],
                                       self.rule['bedtoolstmp']['velcro_bed'])
        if self.debug:
            self.if_runcmd(self.rule['bedtoolstmp'][a_type+"_bed"], cmd)
        else:
            self.run_cmd(cmd)

    def extract(self, a_type):
        """
        extract dhs overlap and velcro overlap information
        """
        self._format(a_type)
        lenall = len(open(self.rule['macsresult']['treat_peaks'], 'rU').readlines())
        if a_type == 'dhs':
            lena_type = len(open(self.rule['bedtoolstmp']['dhs_bed'], 'r').readlines())
        elif a_type == 'velcro':
            lena_type = len(open(self.rule['bedtoolstmp']['velcro_bed'], 'r').readlines())
        self.ratio[a_type] = lena_type
        if lenall != 0:
            self.ratio[a_type + 'percentage'] = float(lena_type)/lenall
        else:
            self.ratio[a_type + 'percentage'] = 0.0001


    def run(self):
        """
        filter macs2 exceptional positions
        then replace the raw output of macs2
        for summits.bed and peaks.bed
        1.use awk '{if ($2 >= 0 && $2 < $3) print}' 6576_peaks.bed > 6576_peaks.bed.temp
          to avoid negative peaks
        2. use bedClip to avoid chromosome out of range use ceas chr len
          bedClip testid_1_peaks.bed /mnt/Storage/data/sync_cistrome_lib/chromLen/hg19.len test.bed
          
        example:
        /opt/bin/wig_correlation_bigwig_input_only.py -d mm9  -s 10  -m mean  --min-score 2  --max-score 50  -r 6576_cor.R 6576_rep1_treat.bw 6576_rep2_treat.bw -l replicate_1 -l replicate_2
        """
        # filter bed files
        def filter(peaks_summits, tmp):
            cmd = "awk '{if ($2 >= 0 && $2 < $3) print}' %s > %s" % (peaks_summits,
                                                                     tmp)
            if self.debug:
                self.if_runcmd(tmp, cmd)
            else:
                self.run_cmd(cmd)
            cmd = "{0} {1} {2} {3}"
            cmd = cmd.format(self.conf['macs']['bedclip'],
                                       tmp,
                                       self.conf['ceas']['chrom_len'],
                                       peaks_summits)
            if self.debug:
                self.if_runcmd(peaks_summits, cmd)
            else:
                self.run_cmd(cmd)

        for rep in range(self.conf['userinfo']['treatnumber']):
            cmd = "awk '{if ($2 >= 0 && $2 < $3) print}' %s > %s"
            
            cmd = cmd % (self.rule['macsresult']['treatrep_peaks'][rep],\
                             self.rule['macstmp']['treatrep_peaks'][rep])
            if self.debug:
                self.if_runcmd(self.rule['macstmp']['treatrep_peaks'][rep], cmd)
            else:
                self.run_cmd(cmd)
            cmd = "{0} {1} {2} {3}"
            cmd = cmd.format(self.conf['macs']['bedclip'],
                             self.rule['macstmp']['treatrep_peaks'][rep],
                             self.conf['ceas']['chrom_len'],
                             self.rule['macsresult']['treatrep_peaks'][rep])
            if self.debug:
                self.if_runcmd(self.rule['macsresult']['treatrep_peaks'][rep], cmd)
            else:
                self.run_cmd(cmd)
        for rep in range(self.conf['userinfo']['treatnumber']):
            cmd = "awk '{if ($2 >= 0 && $2 < $3) print}' %s > %s" % \
                         (self.rule['macsresult']['rep_summits'][rep],
                          self.rule['macstmp']['rep_summits'][rep])
            self.run_cmd(cmd)
            cmd = "{0} {1} {2} {3}"
            cmd = cmd.format(self.conf['macs']['bedclip'],
                             self.rule['macstmp']['rep_summits'][rep],
                             self.conf['ceas']['chrom_len'],
                             self.rule['macsresult']['rep_summits'][rep])
            self.run_cmd(cmd)
        cmd = "awk '{if ($2 >= 0 && $2 < $3) print}' %s > %s" % \
                     (self.rule['macsresult']['treat_peaks'],
                      self.rule['macstmp']['treat_peaks'])
        self.run_cmd(cmd)
        cmd = "{0} {1} {2} {3}"
        cmd = cmd.format(self.conf['macs']['bedclip'],
                                   self.rule['macstmp']['treat_peaks'],
                                   self.conf['ceas']['chrom_len'],
                                   self.rule['macsresult']['treat_peaks'])
        self.run_cmd(cmd)
        cmd = "awk '{if ($2 >= 0 && $2 < $3) print}' %s > %s" % \
                     (self.rule['macsresult']['summits'],
                      self.rule['macstmp']['summits'])
        self.run_cmd(cmd)
        cmd = "{0} {1} {2} {3}"
        cmd = cmd.format(self.conf['macs']['bedclip'],
                                   self.rule['macstmp']['summits'],
                                   self.conf['ceas']['chrom_len'],
                                   self.rule['macsresult']['summits'])
        self.run_cmd(cmd)

        if self.stepcontrol < 3:
            sys.exit(1)
        self.extract('dhs')
        self.extract('velcro') # mouse velcro?
        for k, v in self.ratio.iteritems():
            self.rendercontent['ratios'][k] = v
        self._render()

        if self.conf['userinfo']['treatnumber'] < 2:
            self.log('No replicates, pass the venn diagram and correlation steps')
        else:
            # venn diagram
            if self.conf['userinfo']['treatnumber'] > 3:
                warn('venn diagram support 3 replicates not well')
            cmd = '{0} -t Overlap_of_Replicates {1} {2}'
            cmd = cmd.format(self.conf['venn']['venn_diagram_main'],
                                       ' '.join(self.rule['macsresult']['treatrep_peaks']),
                                       ' '.join(map(lambda x: "-l replicate_" + str(x), 
                                                    xrange(1, self.conf['userinfo']['treatnumber'] + 1)))
                                       )
            self.run_cmd(cmd)
            self.log('venn diagram succeed')
            # correlation plot
            cmd = '{0} -d {1} -s {2} -m {3} --min-score {4} --max-score {5} -r {6} {7} {8}'
            cmd = cmd.format(self.conf['correlation']['wig_correlation_main'],
                                       self.conf['userinfo']['species'],
                                       self.conf['correlation']['wig_correlation_step'],
                                       self.conf['correlation']['wig_correlation_method'],
                                       self.conf['correlation']['wig_correlation_min'],
                                       self.conf['correlation']['wig_correlation_max'],
                                       self.rule['represult']['cor_r'],
                                       ' '.join(self.rule['macsresult']['treatrep_bw']),
                                       ' '.join(map(lambda x: ' -l replicate_' + str(x), xrange(1, self.conf['userinfo']['treatnumber'] + 1))), )
            self.run_cmd(cmd)
            self.log('correlation plot succeed')


class PipeCEAS(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, peaksnumber = 5000, **args):
        """run ceas from top n peaks"""
        super(PipeCEAS, self).__init__(conf, rule, log, **args)
        self.peaks = peaksnumber
        self.stepcontrol = stepcontrol

    def _format(self):
        """
        input : top n peaks bed( filtered by Venn_Cor._format and macs2 bw files
        get peaksnumber ge >= 10 for ceas, or pvalue top n
        # NOTE: generally, ceas peaks number > conservation summits > seqpos summits
        # default top 5000 peaks for ceas
        TODO: filter peaks.xls
        """
        xls = self.rule['macsresult']['peaks_xls']
        bedge = open(self.rule['ceastmp']['ceasge10_bed'], 'w')
        peaksnum = 0
        # fold change
        for line in open(xls, 'rU'):
            if re.search('chr\w', line[0]):
                line = line.strip()
                line = line.split('\t')
                peaksnum += 1
                if peaksnum <= self.peaks and float(line[7]) >= 10:
                    bed_line = line[0] + '\t' + str(int(line[1]) - 1) + '\t' + line[2] + '\t' +\
                        "macs_peaks_" + str(peaksnum) + '\t' + line[8] + '\n'
                    bedge.write(bed_line)
        bedge.close()
        cmd = 'sort -r -g -k 5 %s > %s' % (self.rule['macsresult']['treat_peaks'], self.rule['macstmp']['sortedbed'])
        self.run_cmd(cmd)
        cmd = 'head -n %s %s > %s' % (self.peaks, self.rule['macstmp']['sortedbed'],  self.rule['ceastmp']['ceasp5000'])
        self.run_cmd(cmd)

    def extract(self):
        '''
        merge pdfs
        need to install ghostscript
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=6602_ceas_combined.pdf -f 6602_ceas.pdf 6602_ceas_CI.pdf
        '''
        # cmd = 'gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={0} -f {1} {2}'
        # cmd = cmd.format(self.rule['root']['ceas_pdf'],
        #                            self.rule['ceasresult']['ceaspdf'],
        #                            self.rule['ceasresult']['ceasci']
        #                            )
        # abandon ceas-ex in the future
        
        # cmd = "cp {0} {1}"
        # cmd = cmd.format(self.rule['root']['ceas_pdf'], self.rule['ceasresult']['ceaspdf'])
        # self.run_cmd(cmd)

    def run(self):
        """
        ceasBw and ceas-exon to generate pdf
        and merge together
        added : name option, refgene option, Not added: promotor sizes(--sizes), bipromotor sizes(--bisizes), relative distance(--rel-dist)
        shell examples:
        /usr/local/bin/ceasBW --name 6602_ceas  -g /mnt/Storage/data/RefGene/hg19.refGene  -b peaks_top5000.bed -w 6602_treat.bw -l /mnt/Storage/data/sync_cistrome_lib/chromLen/hg19.len
        /opt/bin/ceas-exon --name 6602_ceas  -g /mnt/Storage/data/RefGene/hg19.refGene  -b peaks_top5000.bed -w 6602_treat.bw
        """
        if self.stepcontrol < 4:
            sys.exit()
        self._format()
        self.log('ceas top n peaks extract done')
        null_unless = lambda trigger, alter :" " if not trigger else alter
        
        sizes_option = null_unless(self.conf['ceas']['ceas_promoter_sizes'],
                                   ' --sizes %s ' % self.conf['ceas']['ceas_promoter_sizes'])
        bisizes_option = null_unless(self.conf['ceas']['ceas_bipromoter_sizes'],
                                     ' --bisizes %s ' % self.conf['ceas']['ceas_bipromoter_sizes'])
        rel_option = null_unless(self.conf['ceas']['ceas_rel_dist'],
                                 ' --rel-dist %s ' % self.conf['ceas']['ceas_rel_dist'])
        gt_option = null_unless(self.conf['ceas']['ceas_genetable_path'],
                                ' -g %s ' % (self.conf['ceas']['ceas_genetable_path']))
        len_option = null_unless(self.conf["ceas"]['chrom_len'],
                                 ' -l %s ' % (self.conf['ceas']['chrom_len']))
        
        cmd = '{0} --name {1} {2} -b {3} -w {4} {5}'
        cmd = cmd.format(self.conf['ceas']["ceas_main"],
                         self.rule['ceastmp']['ceasname'],
                         gt_option + sizes_option + bisizes_option + rel_option,
                         self.rule['ceastmp']['ceasp5000'],
                         self.rule['macsresult']['treat_bw'],
                         len_option)
        if self.debug:
            self.if_runcmd(self.rule['ceastmp']['ceasname']+".pdf", cmd)
        else:
            self.run_cmd(cmd)
        self.extract()
        self.log('ceas succeed')


class PipeConserv(PipeController):
    def __init__(self, conf, rule, a_type, log, stepcontrol, **args):
        """
        use all the peaks from macs2 to draw conservation plot
        """
        super(PipeConserv, self).__init__(conf, rule, log, **args)
        self.type = a_type
        self.stepcontrol = stepcontrol

    def _format(self):
        """
        input: summits bed (filtered by VennCor._format)
        *todo*:get the top n significant peaks for conservation plot
        """
        cmd = 'convert -resize 500x500 -density 50  tmp.pdf {0} | mv {1} {2} '
        cmd = cmd.format(self.rule['conservresult']['conserv_png'],
                         'tmp.R',
                         self.rule['conservresult']['conserv_r'])
        if self.debug:
            self.if_runcmd(self.rule['conservresult']['conserv_r'], cmd)                           
        else:
            self.run_cmd(cmd)

    def extract(self):
        """
        get top n summits from macs2 summits.bed according to ChiLin config
        default 3000 summits
        """
        cmd = 'sort -r -g -k 5 %s | head -n %s > %s ' % (self.rule['macsresult']['summits'], self.conf['conservation']['peaks_num'], self.rule['conservationtmp']['conserv_topn_summits'])
        self.run_cmd(cmd)

    def run(self):
        """
        conservation plot
        to control the region widths for assessing conservation:
            -w 4000 histone, default for TF or Dnase
        shell example:
        /usr/local/bin/conservation_plot.py -t Conservation_at_summits -d /mnt/Storage/data/sync_cistrome_lib/conservation/hg19/placentalMammals -l Peak_summits 6602_summits.bed {-w 4000}
        convert pdf to png
        """
        self.extract()
        if self.stepcontrol < 5:
            sys.exit()
        cmd = '{0} -t Conservation_at_summits -d {1} -l Peak_summits {2} {3}'
        if self.type == 'TF' or self.type == 'Dnase':
            width_option = ' '
        elif self.type == 'Histone':
            width_option = ' -w 4000'   # for histone using width 4000
        else:
            self.log("conservation plot may not support, use default")
        cmd = cmd.format(self.conf['conservation']['conserv_plot_main'],
                                   self.conf['conservation']['conserv_plot_phast_path'],
                                   self.rule['conservationtmp']['conserv_topn_summits'],
                                   width_option)
        if self.debug:
            self.if_runcmd(self.rule['conservresult']['conserv_png'], cmd)
        else:
            self.run_cmd(cmd)
        self._format()
        self.log("conservation plot succeed!")


class PipeMotif(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, **args):
        """pipeline motit part"""
        super(PipeMotif, self).__init__(conf, rule, log, **args)
        self.stepcontrol = stepcontrol

    def _format(self):
        cmd = 'zip -r -q %s results/' % self.rule['motifresult']['seqpos']
        self.run_cmd(cmd)

    def extract(self):
        """
        input: summits.bed
        macs result filtering again:
            generate top n peaks from summits BED file according to p value
            remove chrM from top n summits
            default top 1000 summits.bed, may modify seqpos config
        """
        cmd = 'awk "/^chr[1-22XY]/" %s |sort -r -g -k 5|head -n %s > %s ' % (self.rule['macsresult']['summits'], self.conf['seqpos']['seqpos_top_peaks'], self.rule['motiftmp']['summits_p1000']) 
        self.run_cmd(cmd)

    def run(self):
        """
        get the top 1000 peaks(default)
        shell example:
            /usr/local/bin/MDSeqPos.py -d  -w 600  -p 0.001  -m cistrome.xml  -s hs top1000_summits.bed hg19
        """
        if self.stepcontrol < 6:
            sys.exit()
        self.extract()
        if self.conf['seqpos']['seqpos_width']:
            seqpos_width_option = ' -w  ' +  self.conf['seqpos']['seqpos_width']
        else:
            seqpos_width_option = ' -w 600 '
        if self.conf['seqpos']['seqpos_pvalue_cutoff']:
            seqpos_pvalue_option = ' -p  ' +  self.conf['seqpos']['seqpos_pvalue_cutoff']
        else:
            seqpos_pvalue_option = ' -p 0.001 '
        if self.conf['seqpos']['seqpos_motif_db_selection']:
            seqpos_db_option = ' -m ' +  self.conf['seqpos']['seqpos_motif_db_selection']
        else:
            seqpos_db_option = ' -m transfac.xml,pbm.xml,jaspar.xml,cistrome.xls '
        if self.conf['userinfo']['species'] == 'hg19':
            seqpos_species_option = ' -s hs '
        elif self.conf['userinfo']['species'] == 'mm9':
            seqpos_species_option = ' -s mm '
        elif self.conf['userinfo']['species'] == 'dm4':
            seqpos_species_option = ' -s dm '
        else:
            self.log("MDseqpos not support the species")

        cmd = '{0} -d {1} {2}  {3} {4} {5} {6}'
        cmd = cmd.format(self.conf['seqpos']['seqpos_main'],
                                   seqpos_width_option,
                                   seqpos_pvalue_option,
                                   seqpos_db_option,
                                   seqpos_species_option,
                                   self.rule['motiftmp']['summits_p1000'],
                                   self.conf['userinfo']['species'] )
        if self.debug:
            self.if_runcmd("./results/table_frame.html", cmd)
        else:
            self.run_cmd(cmd)
        self._format()

def package(conf, names, log):
    """
    package all the results in datasetid folder
    """
    bams = glob('*.bam')
    xls = glob('*.xls')
    summits = glob('*_summits.bed')
    peaks = glob('*_peaks.bed')
    bw = glob('*.bw')
    png = glob('*.png')
    cor = glob('*cor*')
    pdf = glob('_ceas_.pdf')
    r = glob('*_ceas_*.R')
    m = glob('*.zip')
    su = glob('*.txt')
    fls = [bams, xls, summits, peaks, bw, png, pdf, r, m, cor, su]
    folder = 'dataset' + conf['userinfo']['datasetid']
    call('mkdir %s' % folder, shell = True)
    for fs in fls:
        for f in fs:
            call('mv %s %s' % (f, folder), shell = True)
    log('package success')
