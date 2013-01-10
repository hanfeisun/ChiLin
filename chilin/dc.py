"""
Data Analysis Workflow
"""
from controller import PipeController
from pkg_resources import resource_filename
from subprocess import call
import os, sys, re

exists = os.path.exists

class PipeGroom(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, **args):
        """ pipeline bam2fastq part"""
        super(PipeGroom, self).__init__(conf, rule, log, **args)
        self.stepcontrol = stepcontrol

    def run(self):
        """use bedtools bamToFastq to replace tophat bamtofastx
        TODO server no such program
        bamToFastq -i x.bam -fq test.fq
        """
        groom_path = lambda x:x.replace(".bam", ".fastq")
        need_groom = lambda x:".bam" in x
        common_cmd = '{0} -i {1} -fq {2}'
        for treat_rep in range(self.conf.basis.treatnumber):
            if need_groom(self.conf.basis.treat[treat_rep]):
                # cmd  = '{0} -q {1} > {2} ' # tophat bam2fastx
                cmd = common_cmd
                cmd = cmd.format(self.conf.bowtie.bam2fq,
                                 self.conf.basis.treat[treat_rep],
                                 groom_path(self.conf.basis.treat[treat_rep]))
                self.log("bam2fastq is processing %s" % (self.conf.basis.treat[treat_rep]))
                self.run_cmd(cmd)
                self.conf.basis.treat[treat_rep] = groom_path(self.conf.basis.treat[treat_rep])
        for control_rep in range(self.conf.basis.controlnumber):
            if need_groom(self.conf.basis.control[control_rep]):
                # cmd  = '{0} -q {1} > {2}'
                cmd = common_cmd
                cmd = cmd.format(self.conf.bowtie.bam2fq,
                                 self.conf.basis.control[control_rep],
                                 groom_path(self.conf.basis.control[control_rep]))
                self.log("bam2fastq is processing %s" % (self.conf.basis.control[control_rep]))
                self.run_cmd(cmd)
                self.conf.basis.control[control_rep] = groom_path(self.conf.basis.control[control_rep])

class PipeBowtie(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, color, **args):
        """pipeline bowtie part"""
        super(PipeBowtie, self).__init__(conf, rule, log, **args)
        self.rendercontent = {}
        self.stepcontrol = stepcontrol
        self.color = color

    def _sam2bam(self, sam, bam):
        """
        using samtools
        convert bowtie sam result to bam
        example: samtools view -bt chrom_len sam bam
        """
        cmd = '{0} view -bt {1} {2} -o {3}'
        cmd = cmd.format(self.conf.samtools.path,
                         self.conf.samtools.chrlen, # chrom_len depends on species
                         sam, bam)
        if self.debug:
            self.ifnot_runcmd(bam, cmd)
        else:
            self.run_cmd(cmd)

    def _extract(self, cnt, files, control = False):
        """
        inputfile : sam files
        cnt : replicates number
        extract bowtie qc information
        sams = [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...] **args
        sams = self.rendercontent
        sams = {'sams': [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...]} *arg
        use awk to speed up a little
        """
        print "working on extracting"
        print cnt
        for sam_rep in range(cnt):
            cmd = ''' time awk -F \'\\t\' -f %s %s > bowtie.tmp ''' \
                  % (resource_filename("chilin", os.path.join("awk", "bowtie_stats.awk")), files[sam_rep])
            if self.debug:
                self.ifnot_runcmd('bowtie.tmp', cmd)
            else:
                self.run_cmd(cmd)
            with open('bowtie.tmp') as f:
                con = f.readlines()
                total_reads = con[0]
                mapped_reads = con[1]
                uniq_read = con[2].rstrip('\n')
                uniq_location = con[3].rstrip('\n')
                usable_percentage = float(con[2])/float(con[0])*100
                print sam_rep
                print self.rule.bowtie.samtreat
                print self.rule.samtools.bamtreat
            if not control:
                self._sam2bam(self.rule.bowtie.samtreat[sam_rep], self.rule.samtools.bamtreat[sam_rep])
            else:
                self._sam2bam(self.rule.bowtie.samcontrol[sam_rep], self.rule.samtools.bamcontrol[sam_rep])
            # uniq_read = len([i for i in reads_dict if reads_dict[i] == 1])
            # uniq_location = len(location_dict)
            # usable_percentage = float(uniq_read)/float(total_reads)*100
            info = { 'name':self.rule.bowtie.samtreat[sam_rep] if not control else self.rule.bowtie.samcontrol[sam_rep],
                     'total': total_reads,
                     'mapped': mapped_reads,
                     'unireads': uniq_read,
                     'uniloc': uniq_location,
                     'percentage' : str(usable_percentage) + '%'}
            self.rendercontent['sams'].append(info)

    def run(self):
        self.rendercontent['sams'] = []
        # color or not, absolid seq
        cmdif = lambda c: '{0} -p {5} -S -C -m {1} {2} {3} {4} ' if c else \
                        '{0} -p {5} -S -m {1} {2} {3} {4} '
        for treat_rep in range(self.conf.basis.treatnumber):
            cmd  = cmdif(self.color)
            cmd = cmd.format(self.conf.bowtie.path,
                             self.conf.bowtie.maxalign,
                             self.conf.bowtie.index,
                             self.conf.basis.treat[treat_rep],
                             self.rule.bowtie.samtreat[treat_rep],
                             self.threads)
            self.log("bowtie is processing %s" % (self.conf.basis.treat[treat_rep]))
            if self.debug:
                self.ifnot_runcmd(self.rule.bowtie.samtreat[treat_rep], cmd)
            else:
                self.run_cmd(cmd)
            self._sam2bam(self.rule.bowtie.samtreat[treat_rep], self.rule.samtools.bamtreat[treat_rep])
        for control_rep in range(self.conf.basis.controlnumber):
            cmd  = cmdif(self.color)
            cmd = cmd.format(self.conf.bowtie.path,
                             self.conf.bowtie.maxalign,
                             self.conf.bowtie.index,
                             self.conf.basis.control[control_rep],
                             self.rule.bowtie.samcontrol[control_rep],
                             self.threads)
            self.log("bowtie is processing %s" % (self.conf.basis.control[control_rep]))
            if self.debug:
                self.ifnot_runcmd(self.rule.bowtie.samcontrol[control_rep], cmd)
            else:
                self.run_cmd(cmd)
            self._sam2bam(self.rule.bowtie.samcontrol[control_rep], self.rule.samtools.bamcontrol[control_rep])
        if exists(self.datasummary) and self.debug:
            print "skip rendering"
            pass
        else:
            self._extract(self.conf.basis.treatnumber, self.rule.bowtie.samtreat, control = False)
            self._extract(self.conf.basis.controlnumber, self.rule.bowtie.samcontrol, control = True)
            self._render()
        self.log("bowtie run successfully")

class PipeMACS2(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, shiftsize, **args):
        """
        MACS step, separately and merge for sorted bam
        shell example:
        macs2 callpeak -B -q 0.01 --keep-dup 1
        --shiftsize=73 --nomodel  -t /Users/testuser/Documents/testchilin/testid_treat_rep2.bam  -c control.bam -n test.bed
        """
        super(PipeMACS2, self).__init__(conf, rule, log, **args)
        # modify as a hidden option in conf
        self.conf = conf
        try:
            if self.conf.macs.model.lower() == 'yes':
                self.model = True
            else:
                self.model = False
        except:
            self.model = False
        self.macsinfo = {}
        self.stepcontrol = stepcontrol
        self.shiftsize = shiftsize
        self.rendercontent = {}

    def _format(self, bdg, tmp, bw):
        """
        filter bdg file and remove over-border coordinates
        convert bdg to bw
        shell example
        /usr/local/bin/intersectBed -a 6523_rep1_treat.bdg -b /opt/bin/chr_limit/chr_limit_mm9.bed -wa -f 1.00 > 6523_rep1_treat.bdg.tm
        /opt/bin/UCSCTools/bedGraphToBigWig 6523_control.bdg.tmp /mnt/Storage/data/Samtool/chromInfo_mm9.txt 6523_control.bw 
        bedGraphTobw chrom len equal to samtools'
        """
        cmd = '{0} intersect -a {1} -b {2} -wa -f 1.00 > {3}'  # bdg filter
        cmd = cmd.format(self.conf.bedtools.path,
                         bdg,
                         self.conf.samtools.chrbed,
                         tmp)
        if self.debug:
            self.ifnot_runcmd(tmp, cmd)
        else:
            self.run_cmd(cmd)
        cmd = '{0} {1} {2} {3} ' # bedGraphTobigwiggle
        cmd = cmd.format(self.conf.macs.bg2bw,
                         tmp,
                         self.conf.samtools.chrlen,
                         bw)
        if self.debug:
            self.ifnot_runcmd(bw, cmd)
        else:
            self.run_cmd(cmd)

    def extract(self):
        """
        get macsinfo of peak calling
        from peaks.xls
        add fc_10 and ratio
        macs2info = {'ratios': {'totalpeak':a, 'total1': 5...}..} *args
        """
        fhd = open(self.rule.macs.peaksxls,"r")
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
        build model
        macs2 callpeak -g hs -B -q 0.01 --keep-dup 1 -t GSM486702_BI.CD34_Primary_Cells.Input.CD34_39661.bed -c GSM486702_BI.CD34_Primary_Cells.Input.CD34_39661.bed  -n test1             
        specify shiftsize
        macs2 callpeak -g hs -B -q 0.01 --keep-dup 1 --nomodel --shiftsize=73 -t GSM486702_BI.CD34_Primary_Cells.Input.CD34_39661.bed -c GSM486702_BI.CD34_Primary_Cells.Input.CD34_39661.bed  -n test1        
        """
        if self.stepcontrol < 2:
            sys.exit(1)
        # Options
        if 'hg' in self.conf.basis.species:
            genome_option = ' -g hs '
        elif 'mm' in self.conf.basis.species:
            genome_option = ' -g mm '
        else:
            genome_option = ' '
        model_option = " " if self.model else " --nomodel --shiftsize=%s " % self.shiftsize
        control_option = " " if self.conf.basis.controlnumber == 0 else ' -c '+ self.rule.samtools.bamcontrolmerge
        # Commands template
        def cmd_merge(input, output):
           if len(input) >1:
               return '{0} merge -f {1}  {2}'.format(self.conf.samtools.path, output, ' '.join(input))
           elif len(input) == 0:
               return ""
           else:
               return 'cp -rf %s %s' % (' '.join(input), output)

        def cmd_callpeak(input_bam, output):
            cmd = '{0} callpeak {1} -B -q 0.01 --keep-dup 1 {2} -t {3} {4} -n {5}'
            return cmd.format(self.conf.macs.path, genome_option,
                              model_option, input_bam,
                              control_option, output)

        # Commands
        cp = lambda orig, dest: self.smart_run(self.cp(orig, dest))
        merge = lambda orig, dest: self.smart_run(cmd_merge(orig, dest), dest)
        callpeak = lambda orig, dest, test: self.smart_run(cmd_callpeak(orig, dest), test)
#        def callpeak(orig, dest, test, exit): # TODO deal with model failure
#            return self.smart_run(cmd_callpeak(orig, dest), test, exit)
        # Start processing
        merge(self.rule.samtools.bamtreat, self.rule.samtools.bamtreatmerge)
        merge(self.rule.samtools.bamcontrol, self.rule.samtools.bamcontrolmerge)

        for repn in range(self.conf.basis.treatnumber):
            callpeak(self.rule.samtools.bamtreat[repn],
                     self.rule.macs.initrep[repn],
                     all([exists(self.rule.macs.initrep[repn]+ '_treat_pileup.bdg'),
                          exists(self.rule.macs.bdgtreatrep[repn])]))

            # convert macs default name to NameRule
            cp(self.rule.macs.initrep[repn]+ '_treat_pileup.bdg',
               self.rule.macs.bdgtreatrep[repn])
            cp(self.rule.macs.initrep[repn] + '_control_lambda.bdg',
               self.rule.macs.bdgcontrolrep[repn])
            self._format(self.rule.macs.bdgtreatrep[repn], 
                         self.rule.macs.bdgtreatreptmp[repn], 
                         self.rule.macs.treatrepbw[repn])
            self._format(self.rule.macs.bdgcontrolrep[repn],
                         self.rule.macs.bdgcontrolreptmp[repn], 
                         self.rule.macs.controlrepbw[repn])
        if self.conf.basis.treatnumber > 1:
            callpeak(self.rule.samtools.bamtreatmerge,
                     self.rule.macs.initmerge,
                     all([exists(self.rule.macs.initmerge + '_treat_pileup.bdg'), \
                          exists(self.rule.macs.bdgtreat)]))

            # convert macs default name to NameRule
            cp(self.rule.macs.initmerge + '_treat_pileup.bdg',
               self.rule.macs.bdgtreat)
            cp(self.rule.macs.initmerge + '_control_lambda.bdg',
               self.rule.macs.bdgcontrol)
            self._format(self.rule.macs.bdgtreat,
                         self.rule.macs.bdgtreattmp,
                         self.rule.macs.treatbw)
            self._format(self.rule.macs.bdgcontrol,
                         self.rule.macs.bdgcontroltmp,
                         self.rule.macs.controlbw)
        else: # Only one replicate, copy the first element of rep list to destination
            origs = [self.rule.macs.peaksrepxls, self.rule.macs.treatreppeaks,
                     self.rule.macs.summitsrep, self.rule.macs.treatrepbw,
                     self.rule.macs.controlrepbw]
            dests = [self.rule.macs.peaksxls, self.rule.macs.treatpeaks,
                    self.rule.macs.summits, self.rule.macs.treatbw, self.rule.macs.controlbw]
            map(lambda orig, dest: cp(orig[0], dest),
                origs, dests)
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
        cmd = '{0} intersect -wa -u -a {1} -b {2} > {3}'
        if a_type == 'dhs':
            cmd = cmd.format(self.conf.bedtools.path,
                             self.rule.macs.treatpeaks,
                             self.conf.bedtools.dhs,
                             self.rule.bedtools.dhs)
            if self.debug:
                self.ifnot_runcmd(self.rule.bedtools.dhs, cmd)
            else:
                self.run_cmd(cmd)
        if a_type == 'velcro':
            cmd = cmd.format(self.conf.bedtools.path,
                             self.rule.macs.treatpeaks,
                             self.conf.bedtools.velcro,
                             self.rule.bedtools.velcro)
            if self.debug:
                self.ifnot_runcmd(self.rule.bedtools.velcro, cmd)
            else:
                self.run_cmd(cmd)

    def extract(self, a_type, species=True):
        """ species indicate test velcro or not
        extract dhs overlap and velcro overlap information
        """
        self._format(a_type)
        lenall = len(open(self.rule.macs.treatpeaks, 'rU').readlines())
        if a_type == 'dhs':
            lena_type = len(open(self.rule.bedtools.dhs, 'r').readlines())
        elif a_type == 'velcro':
            if species:
                lena_type = len(open(self.rule.bedtools.velcro, 'r').readlines())
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
        important modification: -d mm9 => species chrinfo path
        /opt/bin/wig_correlation_bigwig_input_only.py -d mm9_chr_info.path  -s 10  -m mean  --min-score 2  --max-score 50  -r 6576_cor.R 6576_rep1_treat.bw 6576_rep2_treat.bw -l replicate_1 -l replicate_2
        """
        # filter bed files
        def filter(peaks, tmp):
            cmd = "awk '{if ($2 >= 0 && $2 < $3) print}' %s > %s" % (peaks,
                                                                     tmp)
            if self.debug:
                self.ifnot_runcmd(tmp, cmd)
            else:
                self.run_cmd(cmd)
            cmd = "{0} {1} {2} {3}"
            cmd = cmd.format(self.conf.macs.bedclip,
                             tmp,
                             self.conf.ceas.chrlen,
                             peaks)
            if self.debug:
                self.ifnot_runcmd(peaks, cmd)
            else:
                self.run_cmd(cmd)

        for rep in range(self.conf.basis.treatnumber):
            filter(self.rule.macs.treatreppeaks[rep],
                   self.rule.macs.peaksreptmp[rep])
            filter(self.rule.macs.summitsrep[rep],
                   self.rule.macs.summitsreptmp[rep])
        filter(self.rule.macs.treatpeaks,
               self.rule.macs.peakstmp)
        filter(self.rule.macs.summits,
               self.rule.macs.summitstmp)

        if self.stepcontrol < 3:
            sys.exit(1)
        self.extract('dhs', True)
        if self.conf.basis.species == 'hg19':
            self.extract('velcro', True) # mouse has no velcro
        for k, v in self.ratio.iteritems():
            self.rendercontent['ratios'][k] = v
        self._render()
        if self.conf.basis.treatnumber < 2:
            self.log('No replicates, pass the venn diagram and correlation steps')
        else:
            # venn diagram
            if self.conf.basis.treatnumber > 3:
                warn('venn diagram support 3 replicates not well')
            cmd = '{0} -t Overlap_of_Replicates {1} {2} && mv {3} {4}'
            cmd = cmd.format(self.conf.venn.path,
                             ' '.join(self.rule.macs.treatreppeaks),
                             ' '.join(map(lambda x: "-l replicate_" + str(x),
                                          xrange(1, self.conf.basis.treatnumber + 1))),
                             'venn_diagram.png',
                             self.rule.corven.venpng
                            )
            if self.debug:
                self.ifnot_runcmd(self.rule.corven.venpng, cmd)
            else:
                self.run_cmd(cmd)
            self.log('venn diagram succeed')
            # correlation plot
            cmd = '{0} -d {1} -s {2} -m mean --min-score {3} --max-score {4} -r {5} {6} {7} && mv {8}.pdf {9}'
            cmd = cmd.format(self.conf.correlation.path,
                             self.conf.correlation.species,
                             self.conf.correlation.wig_correlation_step,
                             self.conf.correlation.wig_correlation_min,
                             self.conf.correlation.wig_correlation_max,
                             self.rule.corven.corR,
                             ' '.join(self.rule.macs.treatrepbw),
                             ' '.join(map(lambda x: ' -l replicate_' + str(x), xrange(1, self.conf.basis.treatnumber + 1))),
                             self.rule.corven.corR,
                             self.rule.corven.corpdf
                             )
            if self.debug:
                self.ifnot_runcmd(self.rule.corven.corR, cmd)
            else:
                self.run_cmd(cmd)
            self.log('correlation plot succeed')

class PipeCEAS(PipeController):
    def __init__(self, conf, rule, log, stepcontrol, peaksnumber = 5000, **args):
        """run ceas from top n peaks"""
        super(PipeCEAS, self).__init__(conf, rule, log, **args)
        self.peaks = peaksnumber
        self.type = args.get('a_type', 1)
        self.stepcontrol = stepcontrol

    def _format(self):
        """
        input : top n peaks bed( filtered by Venn_Cor._format and macs2 bw files
        get peaksnumber ge >= 10 for ceas, or pvalue top n
        # NOTE: generally, ceas peaks number > conservation summits > seqpos summits
        # default top 5000 peaks for ceas
        TODO: filter peaks.xls
        """
        xls = self.rule.macs.peaksxls
        bedge = open(self.rule.ceas.ge10bed, 'w')
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
        cmd = 'sort -r -g -k 5 %s > %s' % (self.rule.macs.treatpeaks, self.rule.macs.sortedbed)
        self.run_cmd(cmd)

        if self.type == "Dnase": # for Dnase, use all peaks for ceas
            cmd = 'mv %s %s' % (self.rule.macs.sortedbed, self.rule.ceas.peakstop)
        else:   # customized peaks number or default
            cmd = 'head -n %s %s > %s' % (self.peaks, self.rule.macs.sortedbed,  self.rule.ceas.peakstop)
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

        sizes_option = null_unless(self.conf.ceas.promoter_sizes,
                                   ' --sizes %s ' % self.conf.ceas.promoter_sizes)
        bisizes_option = null_unless(self.conf.ceas.bipromoter_sizes,
                                     ' --bisizes %s ' % self.conf.ceas.bipromoter_sizes)
        rel_option = null_unless(self.conf.ceas.rel_dist,
                                 ' --rel-dist %s ' % self.conf.ceas.rel_dist)
        gt_option = null_unless(self.conf.ceas.refgene,
                                ' -g %s ' % (self.conf.ceas.refgene))
        len_option = null_unless(self.conf.ceas.chrlen,
                                 ' -l %s ' % (self.conf.ceas.chrlen))

        cmd = '{0} --name {1} {2} -b {3} -w {4} {5}'
        cmd = cmd.format(self.conf.ceas.path,
                         self.rule.ceas.name,
                         gt_option + sizes_option + bisizes_option + rel_option,
                         self.rule.ceas.peakstop,
                         self.rule.macs.treatbw,
                         len_option)
        if self.debug:
            self.ifnot_runcmd(self.rule.ceas.name+".pdf", cmd)
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
        cmd = 'convert -resize 500x500 -density 50  tmp.pdf {0} && mv {1} {2} '
        cmd = cmd.format(self.rule.conservation.conservpng,
                         'tmp.R',
                         self.rule.conservation.conservR)
        if self.debug:
            self.ifnot_runcmd(self.rule.conservation.conservR, cmd)
        else:
            self.run_cmd(cmd)

    def extract(self):
        """
        get top n summits from macs2 summits.bed according to ChiLin config
        default 3000 summits
        """
        cmd = 'sort -r -g -k 5 %s | head -n %s > %s ' % (self.rule.macs.summits, self.conf.conservation.peaks,
                                                         self.rule.conservation.conservtopsummits)
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
        cmd = cmd.format(self.conf.conservation.path,
                         self.conf.conservation.phast,
                         self.rule.conservation.conservtopsummits,
                         width_option)
        if self.debug:
            self.ifnot_runcmd(self.rule.conservation.conservpng, cmd)
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
        cmd = 'zip -r -q %s results/' % self.rule.motif.seqpos
        self.run_cmd(cmd)

    def extract(self):
        """
        input: summits.bed
        macs result filtering again:
        generate top n peaks from summits BED file according to p value
        remove chrM from top n summits
        default top 1000 summits.bed, may modify seqpos config
        """
        cmd = 'awk "/^chr[1-22XY]/" %s |sort -r -g -k 5|head -n %s > %s ' % \
              (self.rule.macs.summits,
               self.conf.seqpos.summitsnumber, self.rule.motif.summitspeaks1000) 
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
        if self.conf.seqpos.seqpos_width:
            seqpos_width_option = ' -w  ' +  self.conf.seqpos.seqpos_width
        else:
            seqpos_width_option = ' -w 600 '
        if self.conf.seqpos.seqpos_pvalue_cutoff:
            seqpos_pvalue_option = ' -p  ' +  self.conf.seqpos.seqpos_pvalue_cutoff
        else:
            seqpos_pvalue_option = ' -p 0.001 '
        if self.conf.seqpos.db:
            seqpos_db_option = ' -m ' +  self.conf.seqpos.db
        else:
            seqpos_db_option = ' -m transfac.xml,pbm.xml,jaspar.xml,cistrome.xls '
        if self.conf.basis.species == 'hg19':
            seqpos_species_option = ' -s hs '
        elif self.conf.basis.species == 'mm9':
            seqpos_species_option = ' -s mm '
        elif self.conf.basis.species == 'dm4':
            seqpos_species_option = ' -s dm '
        else:
            self.log("MDseqpos not support the species")

        cmd = '{0} -d {1} {2}  {3} {4} {5} {6}'
        cmd = cmd.format(self.conf.seqpos.path,
                         seqpos_width_option,
                         seqpos_pvalue_option,
                         seqpos_db_option,
                         seqpos_species_option,
                         self.rule.motif.summitspeaks1000,
                         self.conf.basis.species)
        if self.debug:
            seqpos_out_path = lambda x:os.path.join("./results",x) # Fixed path
            self.ifnot_runcmd(seqpos_out_path("mdseqpos_out.html"),
                              cmd)
        else:
            self.run_cmd(cmd)
        self._format()

#class PipeReg(PipeController):
#    def __init__(self, conf, rule, log, stepcontrol, **args):
#        """pipeline RegPotential part, for extract the
#        annotated genes near the peaks region"""
#        super(PipeReg, self).__init__(conf, rule, log, **args)
#        self.stepcontrol = stepcontrol
#
#    def _format(self):
#        """
#        REGPOTENTIAL.py -n test_score -t helawithDHS -g /mnt/Storage/data/RefGene/hg19.refGene
#        """
#        pass
#    def run(self):
#        pass

def package(conf, rule, log, **args):
    """
    package all the results in datasetid folder
    """
    from glob import glob
    bams = glob('*.bam')
    xls = glob('*.xls')
    summits = glob('*_summits.bed')
    peaks = glob('*_peaks.bed')
    bw = glob('*.bw')
    png = glob('*.png')
    cor = glob('*cor*')
    pdf = glob('*_ceas.pdf')
    r = glob('*_ceas.R')
    m = glob('*.zip')
    su = glob('dataset*.txt')
    qc = glob('*.tex')
    qcp = glob('*QC.pdf')
    fls = [bams, xls, summits, peaks, bw, png, pdf, r, m, cor, su, qc, qcp]
    folder = 'dataset' + conf.basis.id
    call('mkdir %s' % folder, shell = True)
    for fs in fls:
        for f in fs:
            call('cp %s %s' % (f, folder), shell = True)
    log('package success')
