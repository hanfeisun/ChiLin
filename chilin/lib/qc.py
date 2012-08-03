import os
import math
import re
import zipfile
import subprocess
from pkg_resources import resource_filename
from jinja2 import Environment, FileSystemLoader,PackageLoader


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
    def __init__(self, configs ='',path = ''):
        self.path = path
        self.conf = configs
        self.env = jinja_env
        self.render = {}
        self.summaryCheck = {}
        self.template = self.env.get_template('template.tex')
        self.has_run = False
        self.record = LogWriter().record()


    def run(self):
        """ Run some QC tools or do some time-costing statistics """
        self.has_run = True
        
    
    def _check(self):
        """ Check whether the quality of the dataset is ok. """
        if not self.has_run:
            self.run()
        return True

    def _render(self):
        """ Generate the latex code for current section. """
        temp = self.template.render(self.render)
        self.filehandle.write(temp)
        self.filehandle.flush()


class RawQC(QC_Controller):
    """  
    RawQC aims to perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data.
    """
    def __init__(self,configs = '',path = '', texfile = ''):
        super(RawQC, self).__init__()
        self.conf = configs
        self.path = path
        self.filehandle =  texfile
        
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
        self.has_fastqc = True
        npeakl = []
        nseqlen = []
        for i in range(len(rawdata)):
            d = rawdata[i]
            cmd = '{0} {1} --extract -t 3 -o {2}'
            cmd = cmd.format(self.conf['qc']['fastqc_main'],d,self.path['qcresult']['folder'])
            call(cmd,shell=True)
            temp = os.path.split(d)[1]
            fastqc_out = os.path.join(self.path['qcresult']['folder'],os.path.splitext(temp)[0]+'_fastqc')
            changed_name = os.path.join(self.path['qcresult']['folder'],names[i]+'_fastqc')
            cmd = 'mv {0} {1}'
            cmd = cmd.format(fastqc_out,changed_name)
            call(cmd,shell=True)
            call('rm %s.zip'% fastqc_out,shell=True)
            dataname = changed_name+'/fastqc_data.txt'
            seqlen,peak = self._infile_parse(dataname)
            npeakl.append(peak)
            nseqlen.append(seqlen)
        fastqc_summary = []    #fasqtQC summary
        rCode = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['fastqc_pdf_r'])
        pdfName = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['fastqc_pdf'])
        names = map(lambda x: x.replace('_', '\_'), names)
        for j in range(len(npeakl)):
            if npeakl[j] < 25:
                judge = 'Fail'
            else:
                judge = 'Pass'
            temp = ['%s' % names[j],'%s' % str(nseqlen[j]),'%s' % str(npeakl[j]), '%s' % judge]
            fastqc_summary.append(temp)

        historyData = resource_filename("chilin","db.fastqc_value_list.txt")
        inf=open(historyData,'rU')
        peaklist1=inf.readline()
        f=open(rCode,'w')
        f.write("setwd('%s')\n" %self.conf['userinfo']['outputdirectory'])
        f.write("sequence_quality_score<-c(%s)\n" % str(peaklist1)[1:-1])
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
        inf.close()
        call('Rscript %s' % rCode, shell = True)
        return fastqc_summary, pdfName
    def run(self):
        """ Run some RawQC functions to get final result."""
        self.render['RawQC_check'] = True
        self.render['prefix_datasetid'] = self.conf['userinfo']['datasetid']
        if len(self.conf['userinfo']['controlpath']) ==0:
            rawdata = self.conf['userinfo']['treatpath'].split(',')
            names = self.path['qcresult']['treat_data']
        else:
            rawdata = self.conf['userinfo']['treatpath'].split(',') +self.conf['userinfo']['controlpath'].split(',')
            names = self.path['qcresult']['treat_data']+self.path['qcresult']['control_data']
        for i in range(len(rawdata)-1,-1,-1):
            if '.fastq' in rawdata[i] or '.bam' in rawdata[i] or '.fq' in rawdata[i]:
                pass
            else:
                del rawdata[i]
        if len(rawdata)!=0:
            self.render['fastqc_table'],self.render['fastqc_graph'] = self._fastqc_info(rawdata,names)
            self.render['fastqc_check'] = True
        else:
            self.render['fastqc_check'] = Fasle
        self._check()
        self._render()
    def _check(self):
        """
        Check whether the FastQC's result is ok
        """

         
        
class MappingQC(QC_Controller):
    """ MappingQC aims to describe the mapping quality of the sequence alignment. """
    def __init__(self,configs = '',path = '', texfile = ''):
        super(MappingQC, self).__init__()
        self.conf = configs
        self.path = path
        self.filehandle = texfile

    def _basic_mapping_statistics_info(self,bowtieresult = ''):
        """ Stastic summary of mapping result for each sample. """
        fhd = open(bowtieresult)
        summary = []
        names,totleReads,mappedReads,uniqueReads,uniqueLocation,mapRatio =[],[],[],[],[],[]
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
                mapRatio.append(ratio)
        namesr = map(lambda x: x.replace('_', ' '), names)
        for i in range(len(namesr)):
            temp = [namesr[i],totleReads[i],mappedReads[i],uniqueLocation[i],mapRatio[i]]
            summary.append(temp)
        print summary
        return summary,names,mapRatio

        
    def _mappable_ratio_info(self,ratioList,names):
        """ Cumulative percentage plot to  describe the  mappable ratio quality of all historic data. """
        historyData = self.historyData[0]
        rCode = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['mappable_ratio_r'])
        pdfName = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['mappable_ratio'])
        f=open("%s"% rCode,"w")
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" %pdfName)
        f.write("map_ratio_data<-c(%s)\n" %str(historyData)[0:-1])
        f.write("fn<-ecdf(map_ratio_data)\n")
        f.write("plot(ecdf(map_ratio_data), verticals=TRUE,col.hor='blue', col.vert='black',main='mappable rates',xlab='mappable rates',ylab='Fn(mappable rates)')"+"\n")
        j=0
        for p in ratioList:
            f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(round(p,3),round(p,3),int(pch[j]),col[j]))
            j=j+1
        f.write("legend('topleft',c(%s),pch=c(%s),pt.bg=c(%s))\n" %(str(names)[1:-1],str(pch[:len(names)])[1:-1],str(col[:len(names)])[1:-1]))
        f.write("dev.off()\n")
        f.close()
        cmd = 'Rscript %s'%rCode
        call(cmd,shell=True)
        return pdfName

    def _redundant_ratio_info(self,bamList):
        """ Show redundant  ratio of the dataset in all historic data"""
        names = [os.path.splitext(os.path.split(i)[1])[0] for i in bamList]
        ratioList = []
        for bamfile in bamList:
            temp = 'temp.bed'
            if self.conf['userinfo']['species']=='hg19':
                cmd = 'macs2 filterdup --keep-dup=1 -t {0} -g {1} -o {2}'
                cmd = cmd.format(bamfile,'hs',temp)
                a = subprocess.Popen(cmd,stderr = subprocess.PIPE, shell=True)
                
                content = a.communicate()
                content = list(content)[1].split('\n')
                print content
                os.system('rm temp.bed')
                for line in content:
                    judge = re.findall(r'Redundant rate of alignment file',line)
                    if judge:
                        score = line.split(':')[-1].strip()
                        score = round(1-float(score),3)
                        ratioList.append(score)

        pdfName = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['redundant_ratio'])
        rCode = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['redundant_ratio_r'])
        historyData = self.historyData[2]
        f=open("%s"%rCode,"w")
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" %pdfName)
        f.write("redun_data<-c(%s)\n" % str(historyData)[0:-1])
        f.write("fn<-ecdf(redun_data)\n")
        f.write("plot(ecdf(redun_data), verticals=TRUE,pch='.',main='Unique reads rate ',xlab='Unique reads ratio',ylab='Fn(Unique reads ratio)')"+"\n")
        j=0
        for p in ratioList:
            f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(round(p,3),round(p,3),int(pch[j]),col[j]))
            j=j+1
        f.write("legend('topleft',c(%s),pch=c(%s),pt.bg=c(%s))\n" %(str(names)[1:-1],str(pch[:len(names)])[1:-1],str(col[:len(names)])[1:-1]))
        f.write("dev.off()\n")
        f.close()
        cmd = 'Rscript %s'% rCode
        call(cmd,shell=True)
        if os.path.exists(pdfName):
            return pdfName
        else:
            return Fasle
        
    def run(self,bowtieresult,bampath = ''):
        """ Run some MappingQC function to get final result.
            input: mapping result and path of bam file.  
        """
        self.render['MappingQC_check'] = True
        bampath = os.path.join(self.conf['userinfo']['outputdirectory'],'bowtie')
        bamList = os.popen( "find %s -name \"%s\""%(bampath,'*.bam'))
        bamList = bamList.readlines()
        bamList = [i.strip() for i in bamList]
        historyData = resource_filename("chilin","db.all_data.txt")

        f = open(historyData)
        self.historyData = f.readlines()
        f.close()

        self.render['basic_map_table'],names,mappedRatio = self._basic_mapping_statistics_info(bowtieresult)
        self.render['mappable_ratio_graph'] = self._mappable_ratio_info(mappedRatio,names)
        self.render['redundant_ratio_graph'] = self._redundant_ratio_info(bamList)
        self._render()
        
    def check():
        """Check whether the MappingQC's result is ok. """
        print 'mapping qc pass or not'


class PeakcallingQC(QC_Controller):
    """ PeakcallingQC aims to describe the quality of peak calling result."""
    def __init__(self,configs = '',path = '',texfile = ''):
        super(PeakcallingQC, self).__init__()
        self.conf = configs
        self.path = path
        self.filehandle = texfile
        
    def _peak_summary_info(self,peaksxls):
        """Basic statistic of peak calling result."""
        name = 'dataset'+self.conf['userinfo']['datasetid']
        fhd = open(peaksxls,"rU" )
        float_fc = []
        for i in fhd:
            i = i.strip()
            if i.startswith("# qvalue cutoff"):
                cutoff = i.split('=')[1]
            if i.startswith("# d"):
                shiftsize = int(i.split('=')[1])/2
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
        fhd.close()
        peaks_summary = ['%s'%name,'%s'%cutoff,'%d'%self.totalpeaks,'%d'%self.fold_20,'%s'%shiftsize]
        print peaks_summary
        return peaks_summary
        

    def _high_confidentPeaks_info(self):
        """
        cummulative percentage of peaks foldchange great than 10
        """       
        historyDataName =resource_filename("chilin","db.lg_fold_10.txt")
        pdfName = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['fold_ratio'])
        rCode = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['fold_ratio_r'])
        lg_10 = math.log(self.fold_10,10)
        historyData = open(historyDataName).readlines()[0].strip()
        f = open(rCode,'w')
        f.write('peaks_fc <- c(%s)\n' %historyData)
        f.write('fn <- ecdf(peaks_fc)\n')
        f.write('density <- fn(peaks_fc)\n')
        f.write("pdf('%s')\n" %pdfName)
        f.write("plot(ecdf(peaks_fc),verticals=TRUE,col.hor='blue', col.vert='black',main='Fold 10 peaks Distribution',xlab='lg(number of fold_enrichment>10 peaks)',ylab='Cumulative density function of all public data')\n")
        f.write("points(%f,fn(%f),pch=21,pt.bg=c('#FFB5C5'))\n" % (lg_10,lg_10))
        f.write("legend('topleft',c('ratio of foldchange great than 10'),pch=21,pt.bg=c('#FFB5C5'))\n")
        f.write('abline(v=3,col="red")\n')
        f.write("text(3.5,0,'cutoff=3')\n")
        f.write('dev.off()\n')
        f.close()
        call('Rscript %s' % rCode, shell = True)
        return pdfName
        

    def _velcro_ratio_info(self,peakbed):
        """verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
         The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""
        historyFile = self.conf['venn']['velcro_path']
        overlapped_bed_file = "overlapped_bed_file" # temp
        cmd = '{0} -wa -u -a {1} -b  {2} > {3}'
        cmd = cmd.format(self.conf['bedtools']['intersectbed_main'],
                peakbed,
                historyFile,
                overlapped_bed_file
                )
        call(cmd,shell = True)
        fhd = open(overlapped_bed_file,"r")
        num_overlapped_peaks = len(fhd.readlines())
        fhd.close()

        velcro_ratio = float(num_overlapped_peaks)/self.totalpeaks
        rCode = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['velcro_ratio_r'])
        pdfName = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['velcro_ratio'])
        historyData = self.historyData[3][0:-1]
        f = open(rCode,'w')
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n"% pdfName)
        f.write("rawdata<-c(%s)\n" % historyData)
        f.write("fn<-ecdf(rawdata)\n")
        f.write("plot(ecdf(rawdata), verticals=TRUE,col.hor='blue', col.vert='black',main='velro ratio',xlab='velcro ratio',ylab='Fn(velcro ratio)')\n")
        f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(velcro_ratio,velcro_ratio,int(pch[0]),col[0]))
        f.write("legend('topleft',c('ratio overlap with verlcro'),pch=21,pt.bg=c('#FFB5C5'))\n")
        f.write("dev.off()\n")
        f.close()
        call('Rscript %s' % rCode, shell = True)
        return pdfName
        
    def _DHS_ratio_info(self,peakbed):
        """ DHS ratio indicate the percentage of peaks overlap with DHSs site.
        The function can describe  the particularly dataset's DHS ratio quality of all historic data.
        """
        historyFile = self.conf['venn']['dhs_bed_path']
        overlapped_bed_file = "overlapped_dhs"
        cmd = '{0} -wa -u -a {1} -b  {2} > {3}'
        cmd = cmd.format(self.conf['bedtools']['intersectbed_main'],
                        peakbed,
                        historyFile,
                        overlapped_bed_file
                        )
        call(cmd,shell = True)
        fhd = open(overlapped_bed_file,"r")
        num_overlapped_peaks = len(fhd.readlines())
        dhs_ratio = float(num_overlapped_peaks)/self.totalpeaks
        historyData = self.historyData[2][0:-1]
        rCode = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['dhs_ratio_r'])
        pdfName = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['dhs_ratio'])
        f = open(rCode,'w')
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" % pdfName)
        f.write("rawdata<-c(%s)\n" % historyData)
        f.write("fn<-ecdf(rawdata)\n")
        f.write("plot(ecdf(rawdata), verticals=TRUE,col.hor='blue', col.vert='black',main='overlapped_with_DHSs',xlab='overlapped_with_DHSs',ylab='Fn(overlapped_with_DHSs)')"+"\n")
        f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(dhs_ratio,dhs_ratio,int(pch[0]),col[0]))
        f.write("legend('topleft',c('ratio overlap with DHSs'),pch=21,pt.bg=c('#FFB5C5'))\n")
        f.write("dev.off()\n")
        f.close()
        call('Rscript %s' % rCode, shell = True)
        print num_overlapped_peaks
        print dhs_ratio
        return pdfName
        
    def _replicate_info(self,vennGraph = '',correlationPlot = ''):
        """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""
        self.render['replicte_check'] = True
        self.render['venn_graph'] = vennGraph
        self.render['correlation_graph'] = correlationPlot      

    def run(self,peaksxls,peaksbed):
        """ Run some PeakcallingQC function to get final result. 
            input: peaks bed and excel file.
        """
        self.render['PeakcallingQC_check'] = True
        historyDataName = resource_filename("chilin","db.all_data.txt")

        fph = open(historyDataName)
        self.historyData = fph.readlines()
        fph.close()
        self.render['peak_summary_table'] = self._peak_summary_info(peaksxls)
        self.render['high_confident_peak_graph'] = self._high_confidentPeaks_info()
        self.render['DHS_ratio_graph'] = self._DHS_ratio_info(peaksbed)
        if self.conf['userinfo']['species']=='hg19':
            self.render['verlcro_check'] = True
            self.render['velcro_ratio_graph'] = self._velcro_ratio_info(peaksbed)
        if len(self.conf['userinfo']['treatpath'].split(',')) >= 2:
            vennGraph = os.path.abspath('macs2/'+self.path['represult']['ven_png'])
            correlationPlot = os.path.abspath('macs2/'+self.path['represult']['cor_pdf'])
            self._replicate_info(vennGraph,correlationPlot)
        self._render()
        
    def check():
        """ Check whether PeakcallingQC's result is ok. """
        print 'pass or not'

        
class AnnotationQC(QC_Controller):
    """ AnnotationQC aims to describe the quality of annotations after peak calling. """ 
    def __init__(self,configs = '',path = '', texfile = ''):
        super(AnnotationQC, self).__init__()
        self.conf = configs
        self.path = path
        self.filehandle = texfile
        print 'intialization of function qc'
    def _ceas_info(self,peakxls,ceasCode):
        """ Describe peaks' distribution and relative position. """
        fhd = open( peakxls,"r" )
        list_fc = []
        for i in fhd:
            i = i.strip()
            if i.startswith("# Redundant rate in treatment"):
                temp = i.split(":")
                self.redundant_ratio = str(1-float(temp[1]))
            if i and not i.startswith("#") and not i.startswith("chr\t"):
                fs = i.split("\t")
                fc = fs[7]
                list_fc.append(fc)
        fhd.close()
        
        ceastring = open(ceasCode).read()
        Metaregxcontent = re.findall(r'layout\(matrix\(c\(1, 2, 3, 3, 4, 5\)[^z]*abline\(v=3000\.000000,lty=2,col=c\("black"\)\)', ceastring)[0]
        Pieregxcontent = re.findall(r'# Thus, look at the labels of the pie chart[^z]*# ChIP regions over the genome', ceastring)[0]
        piescript = '\n'.join(Pieregxcontent.split('\n')[4:-3]) + '\n'
        Metascript = '\n'.join(Metaregxcontent.split('\n')[1:]) + '\n'
        # plot 
        rCode = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['ceas_qc_r'])
        Metagene = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['ceas_meta_pdf'])
        Ceasprofile = os.path.join(self.path['qcresult']['folder'],self.path['qcresult']['ceas_profile_pdf'])
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
        f.write("plot(fdd2,type='p',col=2,pch=18,main='Peaks distribution',xlim=c(ma,mi),xlab='Fold change of peaks',ylab='Fn(fold change of peaks)')\n")
        f.write('abline(v=20,lty=2,col=3)\n')
        f.write(piescript)
        f.write('dev.off()\n')
        f.write('# the secend graph \n\n')
        f.write("pdf('%s',height=7,width=5.5)\n" %Ceasprofile)
        f.write("nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), width= c(1,1),height=c(1,1),respect=TRUE)\n")
        f.write(Metascript)
        f.write('dev.off()\n')
        f.close()
        call("Rscript %s"% rCode, shell = True)
        return Metagene,Ceasprofile


    def _DictToList(self,root):
        """extract each node information"""
        result = []
        if "node" not in root.keys():
            return []
        if not root['node']:
            return result
        else:
            result.append(root['node'])
            for each in root['children']:
                result.extend(self._DictToList(each))
            return result

    def _motif_info(self,Zippath):
        """ QC of Sepose. """
        outdir = self.path['qcresult']['folder']
        species = 'Homo sapiens'
        zipFile = zipfile.ZipFile(Zippath)
        data = zipFile.read('results/mdseqpos_out.html') 
        if  not os.path.isdir(outdir+"seqposimg/"):
            cmd = 'unzip ' + Zippath + ' -d ' + outdir + ' \"results/img/*\"'
            call(cmd, shell=True)
            print cmd
        output = []
        content = ""
        for i in data.split('\n'):
            if i.startswith('var mtree'):
                content = i.rstrip().replace('var mtree = ', '')#str type
        exec("mdict = %s" % content)#dict type
        mlist = self._DictToList(mdict)
        for i in mlist:
            if i['zscore'] == 'None':
                i['zscore'] = 65535
                if i['factors'] == []:
                    i['factors'] = ['denovo']
        mlist.sort(key=lambda x:x['zscore'])
        for i in mlist:
            if i['zscore'] < -10 and i['id'].find('observed') > 0 \
                    and species in i['species']:
                        logo = '\includegraphics[angle=0,width=0.28\\textwidth]{%s}'%(outdir+'/results/img/'+i['id']+'_192x120'+'.png')
                        output.append([i['factors'][0].upper(),str(i['hits']),str(i['zscore']),logo])
                        if len(output) >= 7:
                            break
        print output
        print 'motif info\n'
        return output


    def run(self,peaksxls,ceasCode,Zippath,conservatioFile = ''):
        """ Run some AnnotationQC function. """
        self.render['AnnotationQC_check'] =  True
        self.render['meta_gene_graph'],self.render['gene_distribution_graph'] = self._ceas_info(peaksxls,ceasCode)
        if os.path.exists(conservatioFile):
            self.render['conservation_check'] = True
            self.render['conservation_graph'] = conservatioFile
        self.render['motif_table'] = self._motif_info(Zippath)
        self._render()
    def check(self):
        """ Check whether AnnotationQC's result is ok. """
        print 'pass or not'


