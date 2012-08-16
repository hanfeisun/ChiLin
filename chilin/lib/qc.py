import os
import sys
import zipfile
import math
import re
import sqlite3
import subprocess
from subprocess import call
from jinja2 import Environment, FileSystemLoader,PackageLoader
from pkg_resources import resource_filename
from chilin.MotifParser import MotifParser

exists = os.path.exists
notzero = lambda x:os.path.exists(x) and os.path.getsize(x) > 0
def _tospace(x):
    if type(x) == str:
        return x.replace("_"," ")
    return x


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
class QC_Controller(object):
    """
    All the class in the module derives from this class
    """
    def __init__(self, conf, rule, log, texfile, **args):
        self.env = jinja_env
        self.render = {}
        self.checkr = []
        self.summaryCheck = {}
        self.conf = conf
        self.rule = rule
        self.texfile = texfile
        self.log = log
        self.debug = args.get("debug", False)
        self.threads = args.get("threads", 1)  
        self.template = self.env.get_template('template.tex')
        self.db = sqlite3.connect(resource_filename('chilin', 'db/QC_based_infomation.db')).cursor()
        self.shiftsize = args.get("shiftsize", "")


    def run_cmd(self, cmd, exit_ = True, error_handler = lambda :True):
        """
        univeral call shell and judge
        """
        self.log("Run command:\t"+cmd)
        if call(cmd, shell = True):
            # If running command encounters error,
            # call() returns a non-zero value

            result = error_handler()
            if exit_:
                raise
                sys.exit(0)
            else:
                return result
        else:
            return True

    def if_runcmd(self, condition, cmd, else_handler = lambda :True, size_check = True):
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
        
    
    def _check(self):
        """ Check whether the quality of the dataset is ok. """
        if len(self.checkr)!=0:
            for line in self.checkr:
                value = round(float(line[2]),3)
                cutoff = round(float(line[3]),3)
                if value >= cutoff:
                    line.append('pass')
                    self.summarycheck.append(line)
                else:
                    line.append('Fail')
                    self.summarycheck.append(line)

    def _render(self, mode="a"):
        """ Generate the latex code for current section. """
        content = self.template.render(self.render).replace('%','\\%')
        with open(self.texfile, mode) as f:
            f.write(content)



class RawQC(QC_Controller):
    """  
    RawQC aims to perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data.
    """
    def __init__(self,conf = '',rule = '', texfile = '',summarycheck = [],log = '', **args):
        super(RawQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        
    def _infile_parse(self,dataname): # extract information from fastqc result file 
        """ Extract information from fastqc result file. """
        fph = open(dataname)
        data = fph.readlines()
        fph.close()
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
            cmd = cmd.format(self.conf['qc']['fastqc_main'],
                             d,
                             self.conf['userinfo']['outputdirectory'],
                             self.threads)

            
            if self.debug:
                self.if_runcmd(changed_name, cmd)
            else:
               self.run_cmd(cmd)
            cmd = 'cp -rf {0} {1}'
            cmd = cmd.format(fastqc_out,changed_name)
            if self.debug:
                self.if_runcmd(changed_name, cmd)
            else:
                self.run_cmd(cmd)
            self.run_cmd('rm %s.zip'% fastqc_out, exit_=False)
            dataname = changed_name+'/fastqc_data.txt'
            seqlen,peak = self._infile_parse(dataname)
            npeakl.append(peak)
            nseqlen.append(seqlen)
        fastqc_summary = []    #fasqtQC summary
        rCode = self.rule['qcresult']['fastqc_pdf_r']
        pdfName = self.rule['qcresult']['fastqc_pdf']
        names = map(_tospace, names)
        for j in range(len(npeakl)):
            temp = ['%s' % names[j],'%s' % str(nseqlen[j]),'%s' % str(npeakl[j])]
            fastqc_summary.append(temp)
            tempcheck = ['FastqQC','%s' % names[j],'%s' % str(npeakl[j]),25]
            self.checkr.append(tempcheck)


        self.db.execute("select peak_number from fastqc_info_tb")
        fastqc_history = self.db.fetchall()
        historyData = [str(i[0]) for i in fastqc_history]
        historyData = ','.join(historyData)

        f=open(rCode,'w')
        f.write("setwd('%s')\n" %self.conf['userinfo']['outputdirectory'])
        f.write("sequence_quality_score<-c(%s)\n" % historyData)
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("ecdf(sequence_quality_score)->fn\n")
        f.write("fn(sequence_quality_score)->density\n")
        f.write("cbind(sequence_quality_score,density)->fndd\n")
        f.write("fndd2<-fndd[order(fndd[,1]),]\n")
        f.write("pdf('%s')\n" % pdfName)
        f.write("plot(fndd2,type='b',pch=18,col=2,main='Sequence Quality Score Cumulative Percentage',ylab='cummulative density function of all public data')\n")
        j=0
        for p in npeakl:
            f.write("points(%d,fn(%d),pch=%d,bg='%s')\n" %(int(p),int(p),int(pch[j]),col[j]))
            j=j+1
        f.write("legend('topleft',c(%s),pch=c(%s),pt.bg=c(%s))\n" %(str(names)[1:-1],str(pch[:len(names)])[1:-1],str(col[:len(names)])[1:-1]))
        f.write("dev.off()\n")
        f.close()
        
        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        return fastqc_summary, pdfName


    def run(self):
        """ Run some RawQC functions to get final result."""
        self.render['RawQC_check'] = True
        self.render['prefix_datasetid'] = _tospace(self.conf['userinfo']['datasetid'])
        if len(self.conf['userinfo']['controlpath']) ==0:
            rawdata = self.conf['userinfo']['treatpath']
            names = self.rule['qcresult']['treat_data']
        else:
            rawdata = self.conf['userinfo']['treatpath'] +self.conf['userinfo']['controlpath']
            names = self.rule['qcresult']['treat_data']+self.rule['qcresult']['control_data']
        for i in range(len(rawdata)-1,-1,-1):
            if '.fastq' in rawdata[i] or '.bam' in rawdata[i] or '.fq' in rawdata[i]:
                pass
            else:
                del rawdata[i]
        print rawdata
        if len(rawdata)!=0:
            self.render['fastqc_table'],self.render['fastqc_graph'] = self._fastqc_info(rawdata,names)
            self.render['fastqc_check'] = True
        else:
            self.render['fastqc_check'] = Fasle
        self._render("w")
        self._check()

         
        
class MappingQC(QC_Controller):
    """ MappingQC aims to describe the mapping quality of the sequence alignment. """
    def __init__(self,conf = '',rule = '', texfile = '',summarycheck = '',log = '', **args):
        super(MappingQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        self.bowtiesummary = self.rule['root']['data_summary']
        self.log("Preparing mapping QC")

    def _basic_mapping_statistics_info(self,bowtiesummary = ''):
        """ Stastic summary of mapping result for each sample. """
        summary = []
        names,totleReads,mappedReads,uniqueReads,uniqueLocation,mapRatio =[],[],[],[],[],[]
        redundant = []
        mapRatior = []
        with open(bowtiesummary) as fhd:
            for line in fhd:
                line.strip()
                if line.startswith('sam file'):
                    names.append(line.split('=')[1].strip().split('.')[0])
                if line.startswith('total reads'):
                    totle = line.split('=')[1].strip()
                    totleReads.append(totle)
                if line.startswith('mapped reads'):
                    mapped = line.split('=')[1].strip()
                    mappedReads.append(mapped)
                if line.startswith('unique location'):
                    uniqueLocation.append(line.split('=')[1].strip())
                    ratio = round(float(mapped)/float(totle),3)
                    mapRatior.append(str(ratio*100)+'%')
                    mapRatio.append(ratio)

        namesr = map(_tospace, names)
        with open(self.rule['qcresult']['filterdup']) as fredundant:
            for line in fredundant:
                score = round(float(line.strip().split('=')[1]),3)
                redundant.append(score)
                
        for i in range(len(namesr)):
            temp = [namesr[i],totleReads[i],mappedReads[i],mapRatior[i],uniqueLocation[i],redundant[i]]
            tempcheck = ['Unique mappable reads',namesr[i],mappedReads[i],5000000]
            self.checkr.append(tempcheck)
            summary.append(temp)
        print summary
        return summary,namesr,mapRatio

        
    def _mappable_ratio_info(self,ratioList,names):
        """ Cumulative percentage plot to  describe the  mappable ratio quality of all historic data. """
    
        self.db.execute("select map_ratio from mapping_tb")
        mappratio_history = self.db.fetchall()
        historyData = [str(i[0]) for i in mappratio_history]
        historyData = ','.join(historyData)

        rCode = self.rule['qcresult']['mappable_ratio_r']
        pdfName = self.rule['qcresult']['mappable_ratio']
        f=open("%s"% rCode,"w")
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" %pdfName)
        f.write("map_ratio_data<-c(%s)\n" %historyData)
        f.write("fn<-ecdf(map_ratio_data)\n")
        f.write("plot(ecdf(map_ratio_data), verticals=TRUE,col.hor='blue',pch='.',col.vert='black',main='Unique mapped rates',xlab='Unique mapped rates',ylab='Fn(Unique mapped rates)')"+"\n")
        j=0
        for p in ratioList:
            f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(round(p,3),round(p,3),int(pch[j]),col[j]))
            j=j+1
        f.write("legend('topleft',c(%s),pch=c(%s),pt.bg=c(%s))\n" %(str(names)[1:-1],str(pch[:len(names)])[1:-1],str(col[:len(names)])[1:-1]))
        f.write("dev.off()\n")
        f.close()
        cmd = 'Rscript %s'%rCode
        self.run_cmd(cmd)
        return pdfName

    def _redundant_ratio_info(self,bamList):
        """ Show redundant  ratio of the dataset in all historic data"""
        self.log('Processing redundant reads')
        names = [os.path.splitext(os.path.split(i)[1])[0] for i in bamList]
        print names
        ratioList = []
        if notzero(self.rule['qcresult']['filterdup']) and self.debug:    
            self.log("filterdup is skipped because %s exists" % self.rule['qcresult']['filterdup'])
        else:
            with open(self.rule['qcresult']['filterdup'],'w') as fph:
                for i in range(len(bamList)):
                    bamfile = bamList[i]
                    temp = bamfile+".filterdup.temp"
                    temp_out = bamfile +".filterout.temp"
                    cmd = 'macs2 filterdup --keep-dup=1 -t {0} -g {1} -o {2} 2>&1 >/dev/null |tee -a {3}'
                    # print stderr both to screen and file, abandon stdout
                    cmd = cmd.format(bamfile,
                                     self.conf['qc']['filterdup_species'],
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
                            score = round(float(score),3)
                            fph.write('%s=%f\n'%(names[i],score))


                            
        with open(self.rule['qcresult']['filterdup']) as fph:
            for line in fph:
                score = round(float(line.strip().split('=')[1]),3)
                ratioList.append(score)

        pdfName = self.rule['qcresult']['redundant_ratio']
        rCode = self.rule['qcresult']['redundant_ratio_r']

        self.db.execute("select redundant_rate from peak_calling_tb")
        redundant_history = self.db.fetchall()
        historyData = [str(i[0]) for i in redundant_history]
        historyData = [i for i in historyData if i!='null']
        historyData = ','.join(historyData)
        f=open("%s"%rCode,"w")
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" %pdfName)
        f.write("redun_data<-c(%s)\n" % historyData)
        f.write("fn<-ecdf(redun_data)\n")
        f.write("plot(ecdf(redun_data), verticals=TRUE,pch='.',main='Redundant rate ',xlab='Redundant ratio',ylab='Fn(Redundant rate)')"+"\n")
        j=0
        for p in ratioList:
            f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(round(p,3),round(p,3),int(pch[j]),col[j]))
            j=j+1
        f.write("legend('topleft',c(%s),pch=c(%s),pt.bg=c(%s))\n" %(str(names)[1:-1],str(pch[:len(names)])[1:-1],str(col[:len(names)])[1:-1]))
        f.write("dev.off()\n")
        f.close()
        cmd = 'Rscript %s'% rCode
        self.run_cmd(cmd)
        return pdfName
        
    def run(self):
        """ Run some MappingQC function to get final result.
            input: mapping result and path of bam file.  
        """
        print 'mapping qc'
        bowtiesummary = self.bowtiesummary
        historyData = resource_filename('chilin', 'db/all_data.txt')
        f = open(historyData)
        self.historyData = f.readlines()
        f.close()
        self.render['MappingQC_check'] = True
        self.render['Bowtie_check'] = True
        bamList = self.rule['bowtieresult']['bam_treat']+self.rule['bowtieresult']['bam_control']

        self.render['redundant_ratio_graph'] = self._redundant_ratio_info(bamList)
        self.render['basic_map_table'], names, mappedRatio = self._basic_mapping_statistics_info(bowtiesummary)
        print self._basic_mapping_statistics_info(bowtiesummary)
        self.render['mappable_ratio_graph'] = self._mappable_ratio_info(mappedRatio,names)

        self._render()
        self._check()

class PeakcallingQC(QC_Controller):
    """ PeakcallingQC aims to describe the quality of peak calling result."""
    def __init__(self,conf = '',rule = '',texfile = '',summarycheck = '',log = '', **args):
        super(PeakcallingQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        self.peaksxls = self.rule['macsresult']['peaks_xls']
        self.peaksbed = self.rule['macsresult']['treat_peaks']
        self.vennGraph = self.rule['represult']['ven_png']
        self.corrPlot = self.rule['represult']['cor_pdf']
        self.corrR = self.rule['represult']['cor_r']
        
    def _peak_summary_info(self,peaksxls):
        """Basic statistic of peak calling result."""
        name = 'dataset'+self.conf['userinfo']['datasetid']
        with open(peaksxls,"rU" ) as fhd:
            float_fc = []
            for i in fhd:
                i = i.strip()
                cutoff = i.split('=')[1] if i.startswith("# qvalue cutoff")  else "unknown"
                if i and not i.startswith("#") and not i.startswith("chr\t"):
                    fs = i.split("\t")
                    fc = fs[7]
                    float_fc.append(float(fc))
            d = sorted(float_fc)
            d20 = [x for x in d if x >= 20]
            d10 = [x for x in d if x >= 10]
            self.totalpeaks = len(d)+0.01
            self.fold_20 = len(d20)+0.01
            self.fold_10 = len(d10)+0.01
        
        peaks_summary = ['%s'%name,'%s'%cutoff,'%d'%self.totalpeaks,'%d'%self.fold_10,'%s'%self.shiftsize]
        self.checkr.append(['Totle peaks ','%s'%name,'%d'%self.totalpeaks,1000])
        return peaks_summary
        

    def _high_confidentPeaks_info(self):
        """
        cummulative percentage of peaks foldchange great than 10
        """
        name = 'dataset'+self.conf['userinfo']['datasetid']
        self.db.execute("select peak_fc_10 from peak_calling_tb")
        highpeaks_history = self.db.fetchall()
        historyData = [str(math.log(i[0]+0.001,10)) for i in highpeaks_history]
#        historyData = [i for i in historyData if i!='null']
        historyData = ','.join(historyData)

        pdfName = self.rule['qcresult']['fold_ratio']
        rCode = self.rule['qcresult']['fold_ratio_r']
        lg_10 = round(math.log(self.fold_10,10),3)
        pointText = str(lg_10)
        f = open(rCode,'w')
        f.write('peaks_fc <- c(%s)\n' %historyData)
        f.write('fn <- ecdf(peaks_fc)\n')
        f.write('density <- fn(peaks_fc)\n')
        f.write("pdf('%s')\n" %pdfName)
        f.write("plot(ecdf(peaks_fc),verticals=TRUE,col.hor='blue', col.vert='black',main='Fold 10 peaks Distribution',xlab='lg(number of fold_enrichment>10 peaks)',ylab='Cumulative density function of all public data')\n")
        f.write("points(%f,fn(%f),pch=21,bg=c('#FFB5C5'))\n" % (lg_10,lg_10))
        f.write("legend('topleft',c('ratio of foldchange greater than 10 : %s'),pch=21)\n"%pointText)
        f.write('abline(v=3,col="red")\n')
        f.write("text(3.5,0,'cutoff=3')\n")
        f.write('dev.off()\n')
        f.close()
        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        self.checkr.append(['Fold change ','%s'%name,'%d'%lg_10,3])
        return pdfName
        

    def _velcro_ratio_info(self,peakbed):
        """verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
         The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""
        name = 'dataset'+self.conf['userinfo']['datasetid']
        historyFile = self.conf['venn']['velcro_path']
        overlapped_bed_file = "overlapped_bed_file" # temp
        cmd = '{0} -wa -u -a {1} -b  {2} > {3}'
        cmd = cmd.format(self.conf['bedtools']['intersectbed_main'],
                peakbed,
                historyFile,
                overlapped_bed_file
                )
        self.run_cmd(cmd)
        fhd = open(overlapped_bed_file,"r")
        num_overlapped_peaks = len(fhd.readlines())
        fhd.close()

        velcro_ratio = round(float(num_overlapped_peaks)/self.totalpeaks,3)
        rCode = self.rule['qcresult']['velcro_ratio_r']
        pdfName = self.rule['qcresult']['velcro_ratio']

        self.db.execute("select velcro_rate from peak_calling_tb")
        velcro_history = self.db.fetchall()
        historyData = [str(i[0]) for i in velcro_history]
        historyData = [i for i in historyData if i!='null']
        historyData = ','.join(historyData)


        pointText = str(velcro_ratio*100)+'%'
        f = open(rCode,'w')
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n"% pdfName)
        f.write("rawdata<-c(%s)\n" % historyData)
        f.write("fn<-ecdf(rawdata)\n")
        f.write("plot(ecdf(rawdata), verticals=TRUE,col.hor='blue', col.vert='black',main='velro ratio',xlab='velcro ratio',ylab='Fn(velcro ratio)')\n")
        f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(velcro_ratio,velcro_ratio,int(pch[0]),col[0]))
        f.write("legend('topleft',c('ratio overlap with verlcro : %s'),pch=21)\n" %pointText)
        f.write("dev.off()\n")
        f.close()
        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        if velcro_ratio >= 0.1:
            judge = 'pass'
        else:
            judge = 'fail'
        self.summarycheck.append(['Overlap with velcro  ','%s'%name,'%f'%velcro_ratio,0.1,judge])
        return pdfName
        
    def _DHS_ratio_info(self,peakbed):
        """ DHS ratio indicate the percentage of peaks overlap with DHSs site.
        The function can describe  the particularly dataset's DHS ratio quality of all historic data.
        """
        name = 'dataset'+self.conf['userinfo']['datasetid']
        historyFile = self.conf['venn']['dhs_bed_path']
        overlapped_bed_file = "overlapped_dhs"
        cmd = '{0} -wa -u -a {1} -b  {2} > {3}'
        cmd = cmd.format(self.conf['bedtools']['intersectbed_main'],
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
        rCode = self.rule['qcresult']['dhs_ratio_r']
        pdfName = self.rule['qcresult']['dhs_ratio']
        f = open(rCode,'w')
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" % pdfName)
        f.write("rawdata<-c(%s)\n" % historyData)
        f.write("fn<-ecdf(rawdata)\n")
        f.write("plot(ecdf(rawdata), verticals=TRUE,col.hor='blue', col.vert='black',main='overlapped_with_DHSs',xlab='overlapped_with_DHSs',ylab='Fn(overlapped_with_DHSs)')"+"\n")
        f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(dhs_ratio,dhs_ratio,int(pch[0]),col[0]))
        f.write("legend('topleft',c('ratio overlap with DHSs : %s'),pch=21)\n"%pointText)
        f.write("dev.off()\n")
        f.close()
        self.run_cmd('Rscript %s' % rCode, exit_ = False)
        self.checkr.append(['Overlap with DHSs  ','%s'%name,'%f'%dhs_ratio,0.8])
        return pdfName
        
    def _replicate_info(self,vennGraph = '',correlationPlot = '',correlationR = ''):
        """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""
        self.render['replicte_check'] = True
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
            judge = 'pass'
        else:
            judge = 'fail'
        self.summarycheck.append(['Replication QC','%s rep treatment'%m,'%s'%str(cor),'0.6',judge])


    def run(self):
        """ Run some PeakcallingQC function to get final result. 
            input: peaks bed and excel file.
        """
        self.log('Processing PeakcallingQC')
        peaksxls,peaksbed,vennGraph,correlationPlot,correlationR = self.peaksxls,self.peaksbed,self.vennGraph,self.corrPlot,self.corrR
        self.render['PeakcallingQC_check'] = True
        if exists(peaksxls):
            self.render['peak_summary_table'] = map(_tospace, self._peak_summary_info(peaksxls))
            self.render['high_confident_peak_graph'] = self._high_confidentPeaks_info()
        if exists(peaksbed):
            self.render['DHS_ratio_graph'] = self._DHS_ratio_info(peaksbed)
        if self.conf['userinfo']['species']=='hg19':
            self.render['verlcro_check'] = True
            self.render['velcro_ratio_graph'] = self._velcro_ratio_info(peaksbed)
        if len(self.conf['userinfo']['treatpath']) >= 2:
            vennGraph = os.path.abspath('macs2/'+self.rule['represult']['ven_png'])
            correlationPlot = os.path.abspath('macs2/'+self.rule['represult']['cor_pdf'])
            self._replicate_info(vennGraph,correlationPlot,correlationR)
        self._render()
        self._check()
        

        
class AnnotationQC(QC_Controller):
    """ AnnotationQC aims to describe the quality of annotations after peak calling. """ 
    def __init__(self,conf = '',rule = '', texfile = '',summarycheck = '',log = '', **args):
        super(AnnotationQC, self).__init__(conf, rule, log, texfile, **args)
        self.summarycheck = summarycheck
        self.peaksxls = self.rule['macsresult']['peaks_xls']
        self.ceasCode = self.rule['ceasresult']['ceasr']
        self.Zippath = self.rule['motifresult']['seqpos']
        self.conservationFile = self.rule['conservresult']['conserv_png']
        self.conservationR = self.rule['conservresult']['conserv_r']
        print 'intialization of function qc'
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
        
        ceastring = open(ceasCode).read()
        Metaregxcontent = re.findall(r'layout\(matrix\(c\(1, 2, 3, 3, 4, 5\)[^z]*abline\(v=3000\.000000,lty=2,col=c\("black"\)\)', ceastring)[0]
        Pieregxcontent = re.findall(r'# Thus, look at the labels of the pie chart[^z]*# ChIP regions over the genome', ceastring)[0]
        piescript = '\n'.join(Pieregxcontent.split('\n')[4:-3]) + '\n'
        Metascript = '\n'.join(Metaregxcontent.split('\n')[1:]) + '\n'
        # plot 
        rCode = self.rule['qcresult']['ceas_qc_r']
        Metagene = self.rule['qcresult']['ceas_meta_pdf']
        Ceasprofile = self.rule['qcresult']['ceas_profile_pdf']
        f = open(rCode,'w')
        list_fcr = ','.join(list_fc)
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
        f.write("plot(fdd2,type='p',col=2,pch=18,main='Peaks distribution',xlab='Fold change of peaks',ylab='Fn(fold change of peaks)')\n")
        f.write('abline(v=20,lty=2,col=3)\n')
        f.write(piescript)
        f.write('dev.off()\n')
        f.write('# the secend graph \n\n')
        f.write("pdf('%s',height=7,width=5.5)\n" %Ceasprofile)
        f.write("nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), width= c(1,1),height=c(1,1),respect=TRUE)\n")
        f.write(Metascript)
        f.write('dev.off()\n')
        f.close()
        self.run_cmd("Rscript %s"% rCode, exit_ = False)
        return Metagene,Ceasprofile

    def _distance(self,x,y):
        if len(x)!=len(y):
            self.log("warning: x and y has different length")
            x = x[::10][:-1]
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
        fph = open(conservationR)
        for line in fph:
            if re.findall(r'y0<-\S*\)',line):
                value = re.findall(r'y0<-\S*\)',line)[0][6:-1]
                value = value.split(',')
        fph.close()
        value = [float(i) for i in value]
        sumvalue = sum(value)
        value = [i/sumvalue for i in value]
        if atype == 'TF' or atype == 'Dnase':
            histotyDataName = resource_filename("chilin", os.path.join("db", "TFcenters.txt"))
            fph = open(histotyDataName)
            historyData = fph.readlines()
            cutoff = len(historyData)/2
        elif atype == "Histone":
            histotyDataName = resource_filename("chilin", os.path.join("db", "Histone_centers.txt"))
            historyData = fph.readlines()
            fph = open(histotyDataName)
            cutoff = len(historyData)/2
        scoreList = []
        for i in range(len(historyData)):
            temp = historyData[i].strip()
            line = temp.split(' ')

            score = self._distance(value,line)
            scoreList.append(score)
        mindist = round(scoreList.index(min(scoreList)),3)
        if mindist <=cutoff:
            judge = 'pass'
        else:
            judge = 'fail'
        fph.close()
        temp = ['Conservation QC','dataset%s'%self.conf['userinfo']['datasetid'],'%f'%mindist,'K-means cluster','%s'%judge]
        self.summarycheck.append(temp)
        return conservationFile

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

    def get_seqpos(self,Zippath):
        zipFile = zipfile.ZipFile(Zippath)
        data = zipFile.read('results/mdseqpos_out.html')
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
            if i['zscore']<-15 and i['id'].find('observed')>0:
                count += 1
                output.append([('00000%d'%count)[-4:],'|'.join(i['factors']), str(i['zscore']), '|'.join(i['species']), '['+str(i['pssm'])+']',str(i['logoImg']),str(i['hits'])])
        outf = open('seqpose.txt','w')
        outf.write('\t'.join(['id','synonym', 'zscore', 'species', 'pssm','logoImg','hits']) + '\n')
    
        for i in output:
            outf.write('\t'.join(i)+'\n')
        outf.close()
        return 'seqpose.txt'




    def motif_info(self,sqposeTable,Zippath):
        outdir = self.conf['userinfo']['outputdirectory']
        p=MotifParser()
        p.ParserTable(sqposeTable)
        s2 = p.motifs.values()
        logoList = []
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
        for i in s2:
            logo = i['synonym']
            logoList += logo
            denovoNum = logo.count('denovo')
            if denovoNum >=2:
                logo = list(set(logo))
                logo[logo.index('denovo')] = 'denovo::%d'%denovoNum
            cmd = 'unzip ' + Zippath + ' -d ' + outdir + ' \'results/%s\'' %str(i['logoImg'][0])
            logor = os.path.join(outdir,'results/',str(i['logoImg'][0]))
            logorr  = '\includegraphics[angle=0,width=0.28\\textwidth]{%s}'% logor
            tempt = [' '.join(logo),str(i['zscore'][0]),str(i['hits'][0]),logorr]
            os.system(cmd)
            output.append(tempt)
        return output,logoList

    def motif_check(self,logoList,factor):
        factor = factor.upper()
        if factor in logoList:
            judge = 'pass'
        else:
            judge = 'fail'
        temp = ['Motif QC','dataste%s'% self.conf['userinfo']['datasetid'],factor,'-15',judge]
        self.summarycheck.append(temp)


    def run(self,atype):
        """ Run some AnnotationQC function. """
        self.log('#Processing AnnotationQC' )
        peaksxls,ceasCode,Zippath,conservationFile,conservationR = self.peaksxls,self.ceasCode,self.Zippath,self.conservationFile,self.conservationR
        self.render['AnnotationQC_check'] =  True
        if exists(ceasCode):
            print "laila"
            print "NIMEI"
            self.render['ceas_check'] = True
            self.render['meta_gene_graph'],self.render['gene_distribution_graph'] = self._ceas_info(peaksxls,ceasCode)
        if exists(conservationFile) and exists(conservationR):
            self.render['conservation_check'] = True
            self.render['conservation_graph'] = self._conservation_info(conservationR,conservationFile,atype)
        if exists(Zippath):
            tempfile = self.get_seqpos(Zippath)
            motifTable,logoList = self.motif_info(tempfile,Zippath)
            if len(motifTable)>0:
                self.render['motif_check'] = True
                self.render['motif_table'] = motifTable
                self.motif_check(logoList,self.conf['userinfo']['factor'])


        self._render()
        print self.summarycheck


class SummaryQC(QC_Controller):
    """Generate summary report for each QC item and package function"""
    def __init__(self,conf = '',rule = '', log = "", texfile = '', **args):
        super(SummaryQC, self).__init__(conf, rule, log, texfile, **args)
        self.conf = conf
        self.rule = rule

    def run(self,checkList):
        self.render['SummaryQC_check'] = True
        def _prune_id(x):
            if type(x) == str:
                if x.startswith("dataset"):
                    return x[7:]
            return x

        self.render['summary_table'] = map(lambda sub_list: map(lambda x:_prune_id(_tospace(x)),
                                                                sub_list),
                                           checkList)
        print self.render
        self._render()
        cmd = "pdflatex {0}".format(self.texfile)
        self.run_cmd(cmd)
    def packfile(self):
        pass






