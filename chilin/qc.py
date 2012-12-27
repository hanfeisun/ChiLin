import os
import sys
import math
import re
import sqlite3
from subprocess import call
from jinja2 import Environment,PackageLoader
from pkg_resources import resource_filename
from chilin.format import digit_format, percent_format

from chilin.controller import QC_Controller
from chilin.motif.MotifParser import MotifParser
from os.path import exists


has_content = lambda x:os.path.exists(x) and os.path.getsize(x) > 0

def underline_to_space(x):
    if type(x) == str:
        return x.replace("_"," ")
    return x


class RawQC(QC_Controller):
    """  
    RawQC aims to perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data.
    """
    def __init__(self,conf = '',rule = '', texfile = '',summarycheck = [],summaryRender = {},log = '', **args):
        super(RawQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        self.summaryRender = summaryRender
    def _infile_parse(self,dataname): # extract information from fastqc result file 
        """ Extract information from fastqc result file. """
        with open(dataname) as fph:
            data = fph.readlines()
        seqquality = []
        seqlen = 0
        quality_flag = 0
        for line in data:
            if line.startswith('#') or not line:
                continue
            hits = re.findall("Sequence length\t(\d+)", line)
            if hits and len(hits)==1 and not seqlen:
                seqlen = int(hits[0])
            if line.startswith(">>Per sequence quality"):
                quality_flag = 1
            if quality_flag:
                seqquality.append(line)
            if line.startswith(">>END_MODULE") and quality_flag:
                quality_flag = 0
        quality_dict = {}
        for i in seqquality[1:-1]:
            i2 = i.strip().split('\t')
            quality_dict[int(i2[0])] = float(i2[1])
        nt=sum(quality_dict.values())
        n=0
        for item in sorted(quality_dict.items(),key=lambda e:e[0],reverse=True):
            n=n+item[1]
            percentage=n/float(nt)
            if percentage > 0.5:
                peak=int(item[0])
                break
            else:
                continue
        return seqlen,peak
        
        
    def _fastqc_info(self,rawdata,names):
        """ QC analysis of the raw Chip-seq data.
            input: the rawdata list and corresponding lable
        """
        self.log("Begine processing fasctqc")
        self.has_fastqc = True
        npeakl = []
        nseqlen = []
        for i in range(len(rawdata)):
            d = rawdata[i]
            temp = os.path.split(d)[1]
            fastqc_out = os.path.splitext(temp)[0]+'_fastqc'
            changed_name = names[i]
            
            cmd = '{0} {1} --extract -t {3} -o {2}'
            cmd = cmd.format(self.conf.qc.path,
                             d,
                             self.conf.basis.output,
                             self.threads)

            if self.debug:
                self.if_runcmd(not exists(changed_name) and not exists(fastqc_out),
                               cmd)
            else:
               self.run_cmd(cmd)
            cmd = 'cp -rf {0} {1}'
            cmd = cmd.format(fastqc_out,changed_name)
            if self.debug:
                self.if_runcmd(changed_name, cmd)
            else:
                self.run_cmd(cmd, exit_=False)
            self.run_cmd('rm %s.zip'% fastqc_out, exit_=False)
            dataname = changed_name+'/fastqc_data.txt'
            seqlen,peak = self._infile_parse(dataname)
            npeakl.append(peak)
            nseqlen.append(seqlen)
        fastqc_summary = []    #fasqtQC summary
        rCode = self.rule.qc.fastqc_r
        pdfName = self.rule.qc.fastqc_pdf
        names = map(underline_to_space, names)
        for j in range(len(npeakl)):
            temp = ['%s' % names[j],'%s' % str(nseqlen[j]),'%s' % str(npeakl[j])]
            fastqc_summary.append(temp)
            self.checker.append({"desc":"FastqQC",
                                 "data": names[j],
                                 "value":npeakl[j],
                                 "cutoff":25})

        self.db.execute("select peak_number from fastqc_info_tb")
        fastqc_history = self.db.fetchall()
        historyData = [str(i[0]) for i in fastqc_history]
        historyData = ','.join(historyData)
        rnames = str(names)[1:-1]
        npeakl = [str(i) for i in npeakl]
        rnpeakl = ','.join(npeakl)
        print rnpeakl,rnames
        Rrender = {'histroyData':historyData,
                    'value':rnpeakl,
                    'name':rnames,
                    'cutoff':25,
                    'pdfName':pdfName,
                    'main':'Sequence Quality Score Cumulative Percentage',
                    'xlab':'sequence quality score',
                    'ylab':'fn(sequence quality score)',
                    'fastqcCheck':True}

        content = self.Rtemplate.render(Rrender)
        with open(rCode, 'w') as f:
            f.write(content)

        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        return fastqc_summary, pdfName


    def run(self):
        """ Run some RawQC functions to get final result."""
        self.render['RawQC_check'] = True
        self.render['prefix_dataset_id'] = underline_to_space(self.conf.basis.id)
        if len(self.conf.basis.control) ==0:
            rawdata = self.conf.basis.treat[:]
            names = self.rule.qc.treat_data[:]
        else:
            rawdata = self.conf.basis.treat +self.conf.basis.control
            names = self.rule.qc.treat_data + self.rule.qc.control_data

        if len(rawdata)!=0:
            self.render['fastqc_table'],self.render['fastqc_graph'] = self._fastqc_info(rawdata,names)
            self.render['fastqc_check'] = True
        else:
            self.render['fastqc_check'] = False
#        self._render("w")
        self._check()
        self.summaryRender = dict(self.summaryRender.items()+self.render.items())

class MappingQC(QC_Controller):
    """ MappingQC aims to describe the mapping quality of the sequence alignment. """
    def __init__(self,conf = '',rule = '', texfile = '',summarycheck = '',summaryRender = {},log = '', **args):
        super(MappingQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        self.summaryRender = summaryRender
        self.bowtiesummary = self.rule.summary.datasummary
        self.log("Preparing mapping QC")

    def _basic_mapping_statistics_info(self,bowtiesummary = ''):
        """ Stastic summary of mapping result for each sample. """
        summary = []
        names, totalReads,mappedReads,uniqueReads,uniqueLocation,mapRatio =[],[],[],[],[],[]
        redundant = []
        mapRatio = []
        with open(bowtiesummary) as fhd:
            for line in fhd:
                line.strip()
                if line.startswith('sam file'):
                    s = line.split('=')[1].strip().split('.')[0]
                    names.append(s)
                if line.startswith('total reads'):
                    total = line.split('=')[1].strip()
                    totalReads.append(total)
                if line.startswith('mapped reads'):
                    mapped = line.split('=')[1].strip()
                    mappedReads.append(mapped)
                if line.startswith('unique location'):
                    uniqueLocation.append(line.split('=')[1].strip())
                    mapRatio.append(round(float(mapped)/float(total),3))

        with open(self.rule.qc.filterdup) as frddt:
            for line in frddt:
                score = round(float(line.strip().split('=')[1]),3)
                redundant.append(score)

        # formating
        name_sp = map(underline_to_space, names)
        digit_tex = lambda x:digit_format(x,sep="\,") # comma in Tex file should be "/,"
        for i in range(len(name_sp)):
            self.checker.append({"desc": 'Unique mappable reads',
                                 "data": name_sp[i],
                                 "value": mappedReads[i],
                                 "cutoff": 5000000})
            # index problem
            try:
                summary.append([name_sp[i],
                                digit_tex(totalReads[i]),
                                digit_tex(mappedReads[i]),
                                percent_format(mapRatio[i]),
                                digit_tex(uniqueLocation[i]),
                                percent_format(redundant[i])])
            except:
                pass
        for i in range(len(name_sp)):
            self.checker.append({"desc": 'Unique location',
                                 "data": name_sp[i],
                                 "value": uniqueLocation[i],
                                 "cutoff": 5000000})            
            
        return summary,name_sp,mapRatio

        
    def _mappable_ratio_info(self,ratioList,names):
        """ Cumulative percentage plot to  describe the  mappable ratio quality of all historic data. """
    
        self.db.execute("select map_ratio from mapping_tb")
        mappratio_history = self.db.fetchall()
        historyData = [str(i[0]) for i in mappratio_history]
        historyData = ','.join(historyData)
        rCode = self.rule.qc.mappable_ratio_r
        pdfName = self.rule.qc.mappable_ratio
        rnames = str(names)[1:-1]
        ratioList = [str(i) for i in ratioList]
        rratioList = ','.join(ratioList)
        Rrender = {'histroyData':historyData,
                    'value':rratioList,
                    'name':rnames,
                    'cutoff':0.5,
                    'pdfName':pdfName,
                    'main':'Unique mapped rates',
                    'xlab':'Unique mapped rates',
                    'ylab':'fn(Unique mapped rates)',
                    'other':True}
        
        content = self.Rtemplate.render(Rrender).replace('%','\\%')
        with open(rCode, 'w') as f:
            f.write(content)        

        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        return pdfName

    def _redundant_ratio_info(self,bamList):
        """ Show redundant  ratio of the dataset in all historic data"""
        self.log('Processing redundant reads')
        names = [os.path.splitext(os.path.split(i)[1])[0] for i in bamList]
        print names
        ratioList = []# store ratio for plot
        if has_content(self.rule.qc.filterdup) and self.debug:
            self.log("filterdup is skipped because %s exists" % self.rule.qc.filterdup)
        else:
            with open(self.rule.qc.filterdup,'w') as fph:
                for i in range(len(bamList)):
                    bamfile = bamList[i]
                    temp = bamfile+".filterdup.temp"
                    temp_out = bamfile +".filterout.temp"
                    cmd = 'macs2 filterdup --keep-dup=1 -t {0} -g {1} -o {2} 2>&1 >/dev/null |tee -a {3}'
                    #cmd = 'macs2 filterdup --keep-dup=1 -i {0} -g {1} -o {2} 2>&1 >/dev/null |tee -a {3}'
                    # print stderr both to screen and file, abandon stdout
                    cmd = cmd.format(bamfile,
                                     self.conf.qc.species,
                                     temp_out,
                                     temp)
                    if self.debug:
                        self.if_runcmd(temp, cmd)
                    else:
                        self.run_cmd(cmd)
                    with open(temp) as tf:
                        content = tf.readlines()
                    for line in content:
                        judge = re.findall(r'Redundant rate of alignment file',line)
                        if judge:
                            score = line.split(':')[-1].strip()
                            score = 1-round(float(score),3)
                            fph.write('%s=%f\n'%(names[i],score))

        with open(self.rule.qc.filterdup) as fph:
            for line in fph:
                score = round(float(line.strip().split('=')[1]),3)
                ratioList.append(score)
        name_sp = map(underline_to_space, names)
        for i in range(len(name_sp)):
            self.checker.append({"desc": 'Non-Redundant ratio',
                                 "data": name_sp[i],
                                 "value": ratioList[i],
                                 "cutoff": 0.8})

        pdfName = self.rule.qc.redundant_ratio
        rCode = self.rule.qc.redundant_ratio_r

        self.db.execute("select redundant_rate from peak_calling_tb")
        redundant_history = self.db.fetchall()
        historyData = [str(i[0]) for i in redundant_history]
        historyData = [str(1-float(i)) for i in historyData if i!='null']
        historyData = ','.join(historyData)
        rnames = str(names)[1:-1]
        ratioList = [str(i) for i in ratioList]
        rratioList = ','.join(ratioList)
        Rrender = {'histroyData':historyData,
                    'value':rratioList,
                    'name':rnames,
                    'cutoff':0.8,
                    'pdfName':pdfName,
                    'main':'Non-Redundant rate',
                    'xlab':'Non-Redundant rate',
                    'ylab':'fn(Non-Redundant rate)',
                    'other':True}

        content = self.Rtemplate.render(Rrender).replace('%','\\%')
        with open(rCode, 'w') as f:
            f.write(content)        

        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        return pdfName

    def run(self):
        """ Run some MappingQC function to get final result.
            input: mapping result and path of bam file.  
        """
        print 'mapping qc'
        print self.rule.samtools.bamtreat
        if not any(map(lambda x: x.endswith(".bed"), self.rule.samtools.bamtreat)):
            bowtiesummary = self.bowtiesummary
            self.render['MappingQC_check'] = True
            self.render['Bowtie_check'] = True
            bamList = self.rule.samtools.bamtreat + self.rule.samtools.bamcontrol
            self.render['redundant_ratio_graph'] = self._redundant_ratio_info(bamList)

            self.render['basic_map_table'], names,mappedRatio= self._basic_mapping_statistics_info(bowtiesummary)
            self.render['mappable_ratio_graph'] = self._mappable_ratio_info(mappedRatio,names)

#            self._render()
            self._check()
            self.summaryRender = dict(self.summaryRender.items()+self.render.items())

class PeakcallingQC(QC_Controller):
    """ PeakcallingQC aims to describe the quality of peak calling result."""
    def __init__(self,conf = '',rule = '',texfile = '',summarycheck = '',summaryRender = {},log = '', **args):
        super(PeakcallingQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        self.summaryRender = summaryRender
        self.peaksxls = self.rule.macs.peaksxls
        self.peaksbed = self.rule.macs.treatpeaks
        self.vennGraph = self.rule.corven.venpng
        self.corrPlot = self.rule.corven.corpdf
        self.corrR = self.rule.corven.corR

    def _peak_summary_info(self,peaksxls):
        """Basic statistic of peak calling result."""
        name = 'dataset'+self.conf.basis.id
        with open(peaksxls,"rU" ) as fhd:
            float_fc = []
            cutoff = "unknown"
            for i in fhd:
                i = i.strip()
                if i.startswith("# qvalue cutoff"):
                    cutoff = i.split('=')[1] 
                if i and not i.startswith("#") and not i.startswith("chr\t"):
                    fs = i.split("\t")
                    fc = fs[7]
                    float_fc.append(float(fc))
            d = sorted(float_fc)
            d20 = [x for x in d if x >= 20]
            d10 = [x for x in d if x >= 10]
            self.totalpeaks = len(d)+0.01
            self.fold_20 = len(d20)+0.01
            self.fold_10 = len(d10)       
        peaks_summary = ['%s'%name,'%s'%cutoff,'%d'%self.totalpeaks,'%d'%self.fold_10,'%s'%self.shiftsize]
#        self.checker.append({"desc":'Total peaks ',
#                             "data": name,
#                             "value": int(self.totalpeaks),
#                             "cutoff":1000})
        self.checker.append({"desc":'Peaks number with fold change greater than 10X  ',
                             "data": name,
                             "value": int(self.fold_10),
                             "cutoff":1000})

        self.fold_10 = len(d10)+0.01
        return peaks_summary

    def _high_confidentPeaks_info(self):
        """
        cummulative percentage of peaks foldchange great than 10
        """
        name = 'dataset'+self.conf.basis.id
        self.db.execute("select peak_fc_10 from peak_calling_tb")
        highpeaks_history = self.db.fetchall()
#        historyData = [str(math.log(i[0]+0.001,10)) for i in highpeaks_history]
        historyData = [str(math.log(i[0]+0.001,10)) for i in highpeaks_history if i[0] > 0]
#        historyData = [i for i in historyData if i!='null']
        historyData = ','.join(historyData)

        pdfName = self.rule.qc.fold_ratio
        rCode = self.rule.qc.fold_ratio_r
        lg_10 = round(math.log(self.fold_10,10),3)

        Rrender = {'histroyData':historyData,
                    'value':lg_10,
                    'name':"'ratio of fold change upper than 10'",
                    'cutoff':3,
                    'pdfName':pdfName,
                    'main':'High confidence peaks distribution',
                    'xlab':'log(Number of Peaks fold upper than 10)',
                    'ylab':'fn(log(Number of Peaks fold upper than 10))',
                    'other':True}
        
        content = self.Rtemplate.render(Rrender).replace('%','\\%')
        with open(rCode, 'w') as f:
            f.write(content)
        self.run_cmd('Rscript %s' % rCode, exit_ = False)
#        self.checker.append({"desc":'Fold change ',
#                             "data":name,
#                             "value":lg_10,
#                             "cutoff":3})
        return pdfName
        

    def _velcro_ratio_info(self,peakbed):
        """verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
         The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""
        name = 'dataset'+self.conf.basis.id
        historyFile = self.conf.bedtools.velcro
        overlapped_bed_file = "overlapped_bed_file" # temp
        cmd = '{0} intersect -wa -u -a {1} -b  {2} > {3}'
        cmd = cmd.format(self.conf.bedtools.path,
                         peakbed,
                         historyFile,
                         overlapped_bed_file
                         )
        self.run_cmd(cmd)
        fhd = open(overlapped_bed_file,"r")
        num_overlapped_peaks = len(fhd.readlines())
        fhd.close()

        velcro_ratio = round(1-float(num_overlapped_peaks)/self.totalpeaks,3)
        rCode = self.rule.qc.velcro_ratio_r
        pdfName = self.rule.qc.velcro_ratio

        self.db.execute("select velcro_rate from peak_calling_tb")
        velcro_history = self.db.fetchall()
        historyData = [str(i[0]) for i in velcro_history]
        historyData = [str(1-float(i)) for i in historyData if i!='null']
        historyData = ','.join(historyData)


        pointText = str(velcro_ratio*100)+'%'
        Rrender = {'histroyData':historyData,
                    'value':velcro_ratio,
                    'name':"'ratio overlap with non-verlcro:%s'"%pointText ,
                    'cutoff':0.8,
                    'pdfName':pdfName,
                    'main':'non-velcro ratio',
                    'xlab':'non-velcro ratio',
                    'ylab':'fn(non-velcro ratio)',
                    'other':True}

        content = self.Rtemplate.render(Rrender)
        with open(rCode, 'w') as f:
            f.write(content)
        self.run_cmd('Rscript %s' % rCode, exit_ = False)
#        if velcro_ratio >= 0.1:
#            judge = 'Fail'
#        else:
#            judge = 'Pass'
#        self.summarycheck.append(['Overlap with velcro  ','%s'%name,'%f'%velcro_ratio,0.1,judge])
        self.checker.append({"desc":'Overlap with non-velcro ',
                             "data":name,
                             "value":velcro_ratio,
                             "cutoff":0.9})
        return pdfName

    def _DHS_ratio_info(self,peakbed):
        """ DHS ratio indicate the percentage of peaks overlap with DHSs site.
        The function can describe  the particularly dataset's DHS ratio quality of all historic data.
        """
        name = 'dataset'+self.conf.basis.id
        historyFile = self.conf.bedtools.dhs
        overlapped_bed_file = "overlapped_dhs" # tmp
        cmd = '{0} intersect -wa -u -a {1} -b  {2} > {3}'
        cmd = cmd.format(self.conf.bedtools.path,
                        peakbed,
                        historyFile,
                        overlapped_bed_file
                        )
        self.run_cmd(cmd)
        fhd = open(overlapped_bed_file,"r")
        num_overlapped_peaks = len(fhd.readlines())
        dhs_ratio = round(float(num_overlapped_peaks)/self.totalpeaks,3)
        self.db.execute("select DHS_rate from peak_calling_tb ")
        dhs_history = self.db.fetchall()
        historyData = [str(i[0]) for i in dhs_history]
        historyData = [i for i in historyData if i!='null']
        historyData = ','.join(historyData)

        pointText = str(dhs_ratio*100)+'%'
        rCode = self.rule.qc.DHS_ratio_r
        pdfName = self.rule.qc.DHS_ratio

        Rrender = {'histroyData':historyData,
                    'value':dhs_ratio,
                    'name':"'ratio overlap with DHSs:%s'"%pointText ,
                    'cutoff':0.8,
                    'pdfName':pdfName,
                    'main':'Overlapped_with_DHSs',
                    'xlab':'Overlapped_with_DHSs',
                    'ylab':'fn(Overlapped_with_DHSs)',
                    'other':True}

        content = self.Rtemplate.render(Rrender)
        with open(rCode, 'w') as f:
            f.write(content)                

        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        self.checker.append({"desc":'Overlap with DHSs  ',
                             "data":'%s'%name,
                             "value":'%f'%dhs_ratio,
                             "cutoff":0.8})
        return pdfName

    def _replicate_info(self,vennGraph = '',correlationPlot = '',correlationR = ''):
        """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""
        self.render['replicate_check'] = True
        self.render['venn_graph'] = vennGraph
        self.render['correlation_graph'] = correlationPlot
        corrR = self.corrR
        fph = open(correlationR)
        f = open('temCor.r','w')
        m = 0 # count replicate number
        for line in fph:
            if re.match(r"[pc]\d? <- ", line):
                f.write(line)
                m = m+1
        m = m - 2
        for i in range(1,m):
            for j in range(i+1, m+1):
                f.write("print(cor(c[,%s],c[,%s],use='complete.obs'))\n"%(i,j))
        f.close()
        fph.close()
        Rresult = os.popen("Rscript %s"%'temCor.r')
        content = Rresult.readlines()
        Rresult.close()
        os.system('rm temCor.r')
        cors = [round(float(i.split()[1]),3) for i in content]
        if len(cors)==0:
            cor = 0.0001
        else:
            cor = round(sum(cors)/len(cors),3)
        if cor >= 0.6:
            judge = 'Pass'
        else:
            judge = 'Fail'
        self.summarycheck.append(['Replication QC','%s rep treatment'%m,'%s'%str(cor),'0.6',judge])


    def run(self):
        """ Run some PeakcallingQC function to get final result. 
        input: peaks bed and excel file.
        """
        self.log('Processing PeakcallingQC')
        peaksxls,peaksbed,vennGraph,correlationPlot,correlationR = self.peaksxls,self.peaksbed,self.vennGraph,self.corrPlot,self.corrR
        self.render['PeakcallingQC_check'] = True
        if exists(peaksxls):
            self.render['peak_summary_table'] = map(underline_to_space, self._peak_summary_info(peaksxls))
            self.render['high_confident_peak_graph'] = self._high_confidentPeaks_info()
        if exists(peaksbed):
            self.render['DHS_ratio_graph'] = self._DHS_ratio_info(peaksbed)
        if self.conf.basis.species =='hg19':
            self.render['verlcro_check'] = True
            self.render['velcro_ratio_graph'] = self._velcro_ratio_info(peaksbed)
        if self.conf.basis.treatnumber >= 2 :
            vennGraph = os.path.abspath(self.rule.corven.venpng)
            correlationPlot = os.path.abspath(self.rule.corven.corpdf)
            self._replicate_info(vennGraph,correlationPlot,correlationR)
            print vennGraph,correlationPlot
#        self._render()
        self._check()
        self.summaryRender = dict(self.summaryRender.items()+self.render.items())


class AnnotationQC(QC_Controller):
    """ AnnotationQC aims to describe the quality of annotations after peak calling. """ 
    def __init__(self,conf = '',rule = '', texfile = '',summarycheck = '',summaryRender = {},log = '', **args):
        super(AnnotationQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        self.summaryRender = summaryRender
        self.peaksxls = self.rule.macs.peaksxls
        self.ceasCode = self.rule.ceas.R
        self.seqpos_out_path = lambda x:os.path.join("./results",x) # Fixed path
        self.seqpos_stats_out = "seqpos.txt"
        self.conservationFile = self.rule.conservation.conservpng
        self.conservationR = self.rule.conservation.conservR
        print 'initialization of function qc'

    def _ceas_info(self,peakxls,ceasCode):
        """ Describe peaks' distribution and relative position. """
        fhd = open( peakxls,"r" )
        list_fc = []
        with open(peakxls) as fhd:
            for i in fhd:
                i = i.strip()
                if i.startswith("# Redundant rate in treatment"):
                    temp = i.split(":")
                    self.redundant_ratio = str(1-float(temp[1]))
                if i and not i.startswith("#") and not i.startswith("chr\t"):
                    fs = i.split("\t")
                    fc = fs[7]
                    list_fc.append(fc)
        with open(ceasCode) as cf:
            ceastring = cf.read()

        Metaregxcontent = re.findall(r'layout\(matrix\(c\(1, 2, 3, 3, 4, 5\)[^z]*abline\(v=3000\.000000,lty=2,col=c\("black"\)\)', ceastring)[0]
        Pieregxcontent = re.findall(r'# Thus, look at the labels of the pie chart[^z]*# ChIP regions over the genome', ceastring)[0]
        piescript = '\n'.join(Pieregxcontent.split('\n')[4:-3]) + '\n'
        Metascript = '\n'.join(Metaregxcontent.split('\n')[1:]) + '\n'
        # plot 
        rCode = self.rule.qc.ceas_QC_r
        Metagene = self.rule.qc.ceas_meta_pdf
        Ceasprofile = self.rule.qc.ceas_profile_pdf
        list_fcr = ','.join(list_fc)

        with open(rCode,'w') as f:
            f.write("pdf('%s',height=11.5,width=8.5)\n" %Metagene )
            f.write('nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow=TRUE),respect=TRUE)\n')
            f.write('peaks_fc <- c(%s)\n' %list_fcr)
            f.write('fn <- ecdf(peaks_fc)\n')
            f.write('density <- fn(peaks_fc)\n')
            f.write('fdd <- cbind(peaks_fc,density)\n')
            f.write('fdd1 <- fdd[order(fdd[,1],decreasing = TRUE),]\n')
            f.write('fdd2 <- cbind(fdd1[,1],1-fdd1[,2])\n')
            f.write('ma <- max(fdd1[,1])\n')
            f.write('mi <- min(fdd1[,1])\n')
            f.write("plot(fdd2,type='p',col='blue',pch=18,main='Peaks distribution',xlab='Fold change of peaks',ylab='Fn(fold change of peaks)')\n")
#            f.write('abline(v=10,lty=2,col="red")\n')
#            f.write("text(11,0,'cutoff=10')\n")
            f.write(piescript)
            f.write('dev.off()\n')
            f.write('# the secend graph \n\n')
            f.write("pdf('%s',height=7,width=5.5)\n" %Ceasprofile)
            f.write("nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), width= c(1,1),height=c(1,1),respect=TRUE)\n")
            f.write(Metascript)
            f.write('dev.off()\n')

        self.run_cmd("Rscript %s"% rCode, exit_ = False)
        return Metagene,Ceasprofile

    def DictToList(self,root):
        """extract each node information"""
        result = []
        if "node" not in root.keys():
            return []
        if not root['node']:
            return result
        else:
            result.append(root['node'])
            for each in root['children']:
                result.extend(self.DictToList(each))
            return result

    def _distance(self,x,y):
        if len(x)!=len(y):
            print 'error'
        lenght = len(x)
        s=[]
        for i in range(lenght):
            s1=math.pow((float(x[i])-float(y[i])) , 2)
            s.append(s1)
        distance=round(math.sqrt(sum(s)),4)
        return distance

    def _conservation_info(self,conservationR,conservationFile,atype):
        """ For TFcenters data 1,2,3 pass, 4,5,6 fail
            For Histone center data 1,2,3,4 pass, 5,6,7,8 fail.
        """
        with open(conservationR) as fph:
            for line in fph:
                if re.findall(r'y0<-\S*\)',line):
                    value = re.findall(r'y0<-\S*\)',line)[0][6:-1]
                    value = value.split(',')
                elif re.findall(r'x<-c\S*\)',line):
                    xlab = line
            value = [float(i) for i in value]
            sumvalue = sum(value)
            value = [i/sumvalue for i in value]


        if atype == 'TF' or atype == 'Dnase':
            histotyDataName = resource_filename("chilin", os.path.join("db", "TFcenters.txt"))

        elif atype == "Histone":
            histotyDataName = resource_filename("chilin", os.path.join("db", "Histone_centers.txt"))

        with open(histotyDataName) as hf:
                historyData = hf.readlines()
                cutoff = len(historyData)/2

        scoreList = []
        for i in range(len(historyData)):
            temp = historyData[i].strip()
            line = temp.split(' ')

            score = self._distance(value,line)
            scoreList.append(score)
        mindist = scoreList.index(min(scoreList))
        if mindist <=cutoff:
            judge = 'pass'
        else:
            judge = 'fail'
        temp = ['Conservation QC','dataset%s'%self.conf.basis.id, '%f' % round(min(scoreList),3),' ','%s'%judge]
        self.summarycheck.append(temp)
#        print historyData[mindist]
#        print value
        judgevalue = historyData[mindist].strip().split(' ')
        judgevalue = [str(i) for i in judgevalue]
        value = [str(i) for i in value]
        ymax = max(value+judgevalue)
        ymin = min(value+judgevalue)
        rCode = self.rule.qc.conservation_compare_r
        pdfName = self.rule.qc.conservation_compare_pdf
        with open(rCode,'w') as f:
            f.write("pdf('%s',height=8.5,width=8.5)\n" % pdfName)
            f.write("%s\n" % xlab)
            f.write("y1<-c(%s)\n" % ','.join(judgevalue))
            f.write("y2<-c(%s)\n" % ','.join(value))
            f.write("ymax<-(%s)\n" %ymax)
            f.write("ymin<-(%s)\n" %ymin)
            f.write("yquart <- (ymax-ymin)/4\n")
            f.write("plot(x,y2,type='l',col='red',main='Normalized Conservation_at_summits',xlab='Distance from the Center (bp)',ylab='Normalized Average Phastcons',ylim=c(ymin-yquart,ymax+yquart))\n")
            f.write("points(x,y1,type='l',col='blue',lty=2)\n")
            f.write("legend('topleft',c('original group','compared group'),lty=c(1,2),col=c('red', 'blue'))\n")
            f.write("dev.off()\n")
        self.run_cmd("Rscript %s"% rCode, exit_ = False)
        return conservationFile,pdfName


    def get_seqpos(self,cutoff):
        cutoff = float(cutoff)
        with open(self.seqpos_out_path("mdseqpos_out.html")) as mf:
            data = mf.read()
        inf = data.split('\n')
        count = 0
        output = []
        for i in inf:
            if i.startswith('var mtree'):
                data = i.rstrip().replace('var mtree = ','')
        exec('mdict=%s'%data)
        mlist = self.DictToList(mdict)
        for i in mlist:
            if i['zscore'] == 'None':
                i['zscore'] = 65535
            if i['factors'] == []:
                i['factors'] = ['denovo']
        mlist.sort(key=lambda x:x['zscore'])
        for i in mlist:
            if i['zscore'] < cutoff and i['id'].find('observed')>0:
                count += 1
                output.append([('00000%d'%count)[-4:],'|'.join(i['factors']),
                               str(i['zscore']), '|'.join(i['species']),
                               '['+str(i['pssm'])+']',str(i['logoImg']),str(i['hits'])])

        with open(self.seqpos_stats_out, "w") as of:
            of.write('\t'.join(['id','synonym', 'zscore', 'species', 'pssm','logoImg','hits']) + '\n')
            for i in output:
                of.write('\t'.join(i)+'\n')
    def motif_top(self):
        p=MotifParser()
        outdir = self.conf.basis.output
        p.ParserTable(self.seqpos_stats_out)
        s2 = p.motifs.values()
        s2.sort(key=lambda x:x['zscore'][0],reverse=True)
        output = []
        for con in range(5):
            try:
                 i = s2[con]
                 logor = os.path.join(outdir,'results/',str(i['logoImg'][0]))
                 tempt = {'motif':i['synonym'],'hits':str(i['hits'][0]),'Zscore':str(i['zscore'][0]),'logo':logor}
                 output.append(tempt)
            except IndexError:
                pass
        return output

    def motif_info(self,atype):
        outdir = self.conf.basis.output
        p=MotifParser()
        p.ParserTable(self.seqpos_stats_out)
        s2 = p.motifs.values()
        s2.sort(key=lambda x:x['zscore'][0],reverse=True)
        i = 0
        while i<len(s2):
            logo = [s2[i]['synonym'][0]]
            for j in range(len(s2)-1,i,-1):
                if p._Similarity(s2[i]['id'][0],s2[j]['id'][0])[0]>3:
                    logo.append(s2[j]['synonym'][0])
                    id = s2[j]['id'][0]
                    del p.motifs[id]
                    del s2[j]
            s2[i]['synonym'] = logo
            i = i+1
        output = []
        logList = []
        for i in s2:
            logo = i['synonym']
            tt = [logList.append(j) for j in logo]
            denovoNum = logo.count('denovo')
            if denovoNum >=2:
                logo.reverse()
                count = denovoNum
                for lo in range(len(logo)-1,0,-1):
                    if logo[lo]=='denovo':
                        del logo[lo]
                        count -= 1
                        if count==1:
                            break
                logo[logo.index('denovo')] = 'denovo:%d'%denovoNum
                print logo      
            logor = os.path.join(outdir,'results/',str(i['logoImg'][0]))
            tempt = {'motif':logo,'hits':str(i['hits'][0]),'Zscore':str(i['zscore'][0]),'logo':logor}
            output.append(tempt)
            print output
        print '--------------------',logList
        
#        factor = self.conf['meta']['factor'].upper()
#        print factor
#        if atype == 'TF' or atype == 'Dnase':
#            if factor in logList:
#                temp = ['Motif QC','dataset%s'%self.conf['meta']['dataset_id'],factor,' ','pass']
#                self.summarycheck.append(temp)
#            else:
#                temp = ['Motif QC','dataset%s'%self.conf['meta']['dataset_id'],factor,' ','fail']
#                self.summarycheck.append(temp)
        return output


    def run(self,atype):
        """ Run some AnnotationQC function. """
        self.log('#Processing AnnotationQC' )
        peaksxls,ceasCode,conservationFile,conservationR = self.peaksxls,self.ceasCode,self.conservationFile,self.conservationR
        self.render['AnnotationQC_check'] =  True
        if exists(ceasCode):
            self.render['ceas_check'] = True
            self.render['meta_gene_graph'],self.render['gene_distribution_graph'] = self._ceas_info(peaksxls,ceasCode)

        if exists(conservationFile) and exists(conservationR):
            self.render['conservation_check'] = True
#            self.render['conservation_graph'] = conservationFile
            self.render['conservation_graph'],self.render['conservation_compare_graph'] = self._conservation_info(conservationR,conservationFile,atype)

        print self.seqpos_out_path("mdseqpos_out.html")
        print "MAI"
        if exists(self.seqpos_out_path("mdseqpos_out.html")):
            print "MAILA"
            self.get_seqpos(-15)
            motifTable = self.motif_info(atype)
            print motifTable
            if len(motifTable)>0:
                self.render['motif_table'] = motifTable
                self.render['motif_check'] = True
                print '#--------------1'
            else:
                self.get_seqpos(0)
                print '#--------------2'
                motifTable = self.motif_top()
                self.render['motif_table'] = motifTable
                self.render['motif_check'] = True
#        self._render()
        self.summaryRender = dict(self.summaryRender.items()+self.render.items())


class SummaryQC(QC_Controller):
    """Generate summary report for each QC item and package function"""
    def __init__(self,conf = '',rule = '', log = "", texfile = '', **args):
        super(SummaryQC, self).__init__(conf, rule, log, texfile, **args)
        self.conf = conf
        self.rule = rule

    def run(self,checkList,summaryRender,onlyqc):
        self.render['SummaryQC_check'] = True
        self.render['ending'] = True
        self.render = dict(summaryRender.items()+self.render.items())
        print self.render
        def _prune_id(x):
            if type(x) == str:
                if x.startswith("dataset"):
                    return x[7:]
            return x

        self.render['summary_table'] = map(lambda x: map(lambda x:_prune_id(underline_to_space(x)), x), checkList)
        self._render()

        cmd = "pdflatex {0}".format(self.texfile)
        self.run_cmd(cmd)
        self.run_cmd(cmd)       # the pdflatex command should be run twice!
        if onlyqc:
            self.packfile()

    def packfile(self):
        qcfolder = '%s_QCresult' %self.conf.basis.id
        os.system('mkdir %s' %qcfolder)
        for iterm in self.rule.qc:
            print self.rule.qc.iterm
            if isinstance(self.rule.qc.iterm, list):
                for subiterm in self.rule.qc.iterm:
                    if has_content(subiterm):
                        cmd = 'cp -rf %s %s'%(subiterm, qcfolder)
                        self.run_cmd(cmd)
                        print subiterm
            elif has_content(self.rule.qc.iterm):
                print self.rule.qc.iterm
                cmd = 'cp -rf %s %s'%(self.rule.qc.iterm, qcfolder)
                self.run_cmd(cmd)



