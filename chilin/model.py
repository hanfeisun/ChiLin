"""
mapping conf and rule dictionaries into simple and
universe programming names method
"""
from jinja2 import Environment, PackageLoader
from meta import Section, RuleBase, ConfBase
from pkg_resources import resource_filename

import logging, os
import sys



def gen_conf( species, platform ):
    """ options for subparser gen """
    env = Environment(loader = PackageLoader('chilin', 'conf'))
    if platform == "Illumina":
        new_temp = "ChiLin.conf"
    elif platform == "absolid":
        new_temp = "ChiLin_ab.conf"
    temp = env.get_template(new_temp)
    if species == 'hg19':
        conf = temp.render(species = species,
                           filterdup_species = 'hs')
    elif species == 'mm9':
        conf = temp.render(species = species,
                           filterdup_species = 'mm')
    print conf

class LogWriter:
    def __init__(self, logfile = 'log'):
        """ Universal log format time + incidence """
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

class Conf(ConfBase):
    basis = Section(["user", "time", "id", "species", "factor", "treat", "control", "output"])
    bowtie = Section(["path", "index", "bam2fq", "maxalign"])
    samtools = Section(["path", "chrlen", "chrbed"])
    macs = Section(["path", "bg2bw", "bedclip"])
    bedtools = Section(["path", "dhs", "velcro"])
    ceas = Section(["path", "exon", "peaks", "refgene", "chrlen", "promoter_sizes", "bipromoter_sizes", "rel_dist"])
    conservation = Section(["path", "peaks", "width", "phast"])
    correlation = Section(["path", "species", "wig_correlation_step", "wig_correlation_min", "wig_correlation_max"])
    venn = Section(["path"])
    seqpos = Section(["path", "species", "summitsnumber", "mdscan_width", "mdscan_top_peaks",
                      "seqpos_mdscan_top_peaks_refine", "seqpos_width", "seqpos_pvalue_cutoff", "db"])
    qc = Section(["path", "species"])

class Rule(RuleBase):
    bowtie = Section(["samtreat", "samcontrol"])
    samtools = Section(["bamtreat", "bamcontrol", "bamtreatmerge", "bamcontrolmerge"])
    bedtools = Section(["dhs", "velcro", "overlap"])
    macs = Section(["initrep", "initmerge",
                    # tmp
                    "sortedbed",
                    "peaksreptmp", "peakstmp",
                    "summitsreptmp", "summitstmp",
                    "bdgcontrolrep", "bdgcontrolreptmp","bdgcontrol","bdgcontroltmp",
                    "bdgtreatrep", "bdgtreatreptmp", "bdgtreat", "bdgtreattmp",
                    "treatreppeaks", "treatpeaks", "treatrepbw", "treatbw",
                    "controlrepbw", "controlbw",
                    "peaksrepxls", "peaksxls",
                    "summitsrep", "summits"
                    ])
    ceas = Section(["name",
                    # tmp
                    # automatically generated tmp
                    # "withoutpeakR", "withpeakpdf",
                    # "withpeakR", "withpeakpdf",
                    "ge10bed", "peakstop",
                    # result
                    # automatically generated results
                    # "CIR", "CIpdf",
                    "R"
                    # "xls", "R", "pdf",
                    # "combined"
                    ])
    corven = Section(["corR", "corpdf",
                      "venR", "venpng"])
    conservation = Section(["conservtopsummits",
                            "conservR", "conservpng"])
    motif = Section(["summitspeaks1000", "seqpos"])
    qc = Section(["treat_data", "control_data", "fastqc_r", "fastqc_pdf", #tmp
                  "filterdup",
                  "mappable_ratio_r", "mappable_ratio",
                  "redundant_ratio_r", "redundant_ratio",
                  "velcro_ratio_r", "velcro_ratio",
                  "fold_ratio_r", "fold_ratio",
                  "DHS_ratio_r", "DHS_ratio",
                  "ceas_QC_r", "ceas_meta_pdf", "ceas_profile_pdf",
                  "conservation_compare_r", "conservation_compare_pdf",
                  # result
                  "QCtex", "QCreport"])
    summary = Section(['datasummary', 'log'])

class PipePreparation(Conf, Rule):
    def __init__(self, ChiLinconfPath,
                 NameConfPath = resource_filename("chilin", os.path.join("conf", "NameRule.conf"))):
        """
        Parse the Name Rule and
        Customer filled Chilinconf
        """
        self.ChiLinconfPath = ChiLinconfPath
        self.NameConfPath = NameConfPath
        self.checked = False
        self._conf = Conf(ChiLinconfPath)
        parseinput = lambda x: [i for i in x.strip().split(',') if i]
        # need absolute path
        os.chdir(self._conf.basis.output)
        self._conf.basis.treat = parseinput(self._conf.basis.treat)
        self._conf.basis.control = parseinput(self._conf.basis.control)
        self._conf.basis.treatnumber = len(self._conf.basis.treat)
        self._conf.basis.controlnumber = len(self._conf.basis.control)
        print self._conf.basis.treatnumber
        print self._conf.basis.treat
        print self._conf.basis.control
        print self._conf.basis.controlnumber
        ## don't disorder control and treat !!
        self._rule = Rule(NameConfPath, self._conf.basis.id,
                          self._conf.basis.controlnumber, self._conf.basis.treatnumber)
        self.log = LogWriter(self._rule.summary.log).record

    def get_config(self):
        if self.checked == False:
            self.check_conf()
        return self._conf

    def get_rule(self):
        if self.checked == False:
            self.check_conf()
        return self._rule

    def get_tex(self):
        return self._rule.qc.QCtex

    def get_summary(self):
        return self._rule.summary.datasummary

    def check_conf(self):
        """
        Check the Meta configuration
        if up to our definition
        """
        is_null = lambda a_list: len(a_list) > 0
        exists = os.path.exists
        error = logging.error
        warn = logging.warning

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
        c.append(fatal(is_null(self._conf.basis.treat),
                       'No treat file path'))
        c.append(danger(is_null(self._conf.basis.control),
                        'No control file path'))
        c.append(fatal(is_null(self.basis.id),
                       'No datasetID'))
        c.append(fatal(is_null(self.basis.factor),
                       'No factor'))
        c.append(fatal(not os.path.isdir(self._conf.basis.output),
                       'check output directory name'))
        c.append(fatal(self._conf.basis.species not in ['hg19', 'mm9'],
                       'No species or current species %s not supported' % self._conf.basis.species))
        c.append(fatal(all(map(os.path.isfile, self._conf.basis.treat)),
                       'check your treat file, one of them is not a file'))
        c.append(fatal(all(map(os.path.isfile, self._conf.basis.control)),
                       'check your control file, one of them is not a file'))
        c.append(fatal(not exists(self._conf.qc.path),
                       'fastqc not exists'))
        c.append(fatal(any(map(lambda x:x.endswith(".fastq") or x.endswith(".bam"),
                               self._conf.basis.treat)),
                       'check your treat file, one of them is not ended with .fastq or .bam'))
        c.append(fatal(any(map(lambda x:x.endswith(".fastq") or x.endswith(".bam") or x.endswith(".bed"),
                               self._conf.basis.control)),
                       'check your control file, one of them is not ended with .fastq or .bam'))                               
        c.append(fatal(not exists(self._conf.bowtie.path),
                       "bowtie program dependency has problem"))
        c.append(fatal(not exists(self._conf.samtools.path),
                       "samtools dependency has problem"))
        c.append(fatal(not exists(self._conf.macs.path),
                       "macs2 program dependency has probelm"))
        c.append(fatal(not exists(self._conf.bedtools.path),
                       "bedtools program dependency has problem"))
        c.append(fatal(not exists(self._conf.conservation.path),
                       "conservation_plot dependency has problem"))
        c.append(fatal(not exists(self._conf.correlation.path),
                       "correlation plot dependency has problem"))
        c.append(fatal(not exists(self._conf.venn.path),
                       "venn plot dependency has problem"))
        check_list(c)
        self.checked = True
        self.log('ChiLin config parse success, and dependency check has passed!')
