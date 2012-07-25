import chilin
from chilin.dc import *
import os
import math
import re
import zipfile
from jinja2 import Environment, FileSystemLoader,PackageLoader
import ChiLin
# FileSystemLoader('/Users/Samleo/mybin/chilin/chilin/lib/template/')
# JinJa_temp = os.path.dirname(os.path.abspath(chilin.__file__))
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
        self.template = self.env.get_template('template.tex')
        self.has_run = False
        print 'QC control'

    def run(self):
        """ Run some QC tools or do some time-costing statistics """
        self.has_run = True
        return True
        
    def _partion(self):
        pass
    
    def _check(self):
        """ Check whether the quality of the dataset is ok. """
        if not self.has_run:
            self.run()
        return True

    def _render(self):
        """ Generate the latex code for current section. """
#       self.template.render({})
        pass


class RawQC(QC_Controller):
    """  
    RawQC aims to perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data.
    """
    def __init__(self,configs = '',path = '', texfile = ''):
#        self.env = jinja_env
#       self.template = self.env.get_template('template.tex')
        super(RawQC, self).__init__()
        self.conf = configs
        self.path = path
        self.filehandle =  texfile
        
    def _infile_parse(self,dataname): # extract information from fastqc result file
        data = open(dataname).readlines()
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
        data.close()
        
    def _fastqc_info(self,rawdata,names):
        self.has_fastqc = True
        """ QC analysis of the raw Chip-seq data, including sequence quality score of particularity raw data and the cumulative percentage plot of the sequence quality scores of all historic data.
        """
        npeakl = []
        nseqlen = []
        for i in range(len(rawdata)):
            d = rawdata[i]
            cmd = '{0} {1} --extract -t 3 -o {2}'
            cmd = cmd.format(self.conf['qc']['fastqc_main'],d,self.path['qcresult']['folder'])
            call(cmd,shell=True)
            tem = d.split('/')[-1]
            fastqc_out = self.path['qcresult']['folder']+'/'+tem.split('.')[0]+'_fastqc'
            changed_name = self.path['qcresult']['folder']+'/'+names[i]+'_fastqc'
            cmd = 'mv {0} {1}'
            cmd = cmd.format(fastqc_out,changed_name)
            call(cmd,shell=True)
            call('rm %s.zip'% fastqc_out,shell=True)
            dataname = changed_name+'/fastqc_data.txt'
            print dataname
            seqlen,peak = self._infile_parse(dataname)
            npeakl.append(peak)
            nseqlen.append(seqlen)
        fastqc_summary = []    #fasqtQC summary
        rCode = self.path['qcresult']['folder']+'/'+self.path['qcresult']['fastqc_pdf_r']
        pdfName = self.path['qcresult']['folder']+'/'+self.path['qcresult']['fastqc_pdf']
        for j in range(len(npeakl)):
            if npeakl[j] < 25:
                judge = 'Fail'
            else:
                judge = 'Pass'
            temp = '%s: %s/t%s/t%s ' %(names[j],str(nseqlen[j]),str(npeakl[j]),judge)
            fastqc_summary.append(temp)
        print fastqc_summary
        historyData = os.path.split(chilin.__file__)[0] + '/' + 'db/fastqc_value_list.txt'
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
        self.RawQC_check = True
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
            self.fastqc_summary_stat,self.fastqc_graph_stat = self._fastqc_info(rawdata,names)
            self.fastqc_check = True
        else:
            self.fastqc_check = Fasle
        self._check()
        self._render()
        return 'test'
    def _check(self):
        """
        Check whether the FastQC's result is ok
        """
        if self.has_fastqc:
            print 'input self.fastqc_summary_stat for judge'
        else:
            print ' no need fastqc'
    def _render(self):
        print self.fastqc_summary_stat
        temp = self.template.render(RawQC_check = self.RawQC_check,prefix_datasetid = 'id',fastqc_check = self.has_fastqc,fastqc_graph = self.fastqc_graph_stat)
        self.filehandle.write(temp)
        self.filehandle.flush()
         


        
class MappingQC(QC_Controller):
    """ MappingQC aims to describe the mapping quality of the sequence alignment. """
    def __init__(self,configs = '',path = '', texfile = ''):
        super(MappingQC, self).__init__()
        self.conf = configs
        self.path = path
        self.filehandle = texfile
        print 'init mapping qc'

    def _basic_mapping_statistics_info(self,mapDict = ''):
        """ Stastic summary of mapping result for each sample. """
        self.mappable_summary_stat = []
        temp =  'this part i still under construct'
        self.mappable_summary_stat.append(temp)
        print self.mappable_summary_stat.append
        
        """ Cumulative percentage plot to  describe the  mappable ratio quality of all historic data. """
    def _mappable_ratio_info(self,ratioList,names):
        historyData = self.historyData[0]
        rCode = self.path['qcresult']['folder']+'/'+self.path['qcresult']['mappable_ratio_r']
        pdfName = self.path['qcresult']['folder']+'/'+self.path['qcresult']['mappable_ratio']
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

        """ Cumulative percentage plot to  describe the  mappable ratio quality of all historic data."""
        print 'mappable_ratio'

    def _redundant_ratio_info(self,ratioList,names):
        """ Show redundant  ratio of the dataset in all historic data"""
        print 'redundant_ratio\n'
        pdfName = self.path['qcresult']['folder']+'/'+self.path['qcresult']['redundant_ratio']
        rCode = self.path['qcresult']['folder']+'/'+self.path['qcresult']['redundant_ratio_r']
        historyData = self.historyData[2]
        f=open("%s"%rCode,"w")
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" %pdfName)
        f.write("redun_data<-c(%s)\n" % str(historyData)[0:-1])
        f.write("fn<-ecdf(redun_data)\n")
        f.write("plot(ecdf(redun_data), verticals=TRUE,pch='.',main='redundant ratio',xlab='redundant ratio',ylab='Fn(redundant ratio)')"+"\n")
        j=0
        for p in ratioList:
            f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(round(p,3),round(p,3),int(pch[j]),col[j]))
            j=j+1
        f.write("legend('topleft',c(%s),pch=c(%s),pt.bg=c(%s))\n" %(str(names)[1:-1],str(pch[:len(names)])[1:-1],str(col[:len(names)])[1:-1]))
        f.write("dev.off()\n")
        f.close()
        cmd = 'Rscript %s'% rCode
        call(cmd,shell=True)
        return pdfName

    def _render(self):
        temp = self.template.render(MappingQC_check = self.MappingQC_check, basic_mapping_table = self.mappable_summary_stat,\
                mappable_ratio_graph = self.mappable_ratio_stat)
        self.filehandle.write(temp)
        self.filehandle.flush()
    def run(self,mapDict):
        
        self.MappingQC_check = True
        names = ['rep1']
        mapped_ratio = [0.8]
        historyData = os.path.split(chilin.__file__)[0] + '/' + 'db/all_data.txt'
        f = open(historyData)
        self.historyData = f.readlines()
        f.close()
        """ Run some MappingQC function to get final result. """
        self._basic_mapping_statistics_info()
        self.mappable_ratio_stat = self._mappable_ratio_info(mapped_ratio,names)
#        self.redundant_ratio_stat = self._redundant_ratio_info(redundant_ratio,names)
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
        print 'init peak calling  qc'
        
    def _peak_summary_info(self,peaksxls,fold = 20):
        """Basic statistic of peak calling result."""
        print ' summary_info'
        fhd = open(peaksxls,"rU" )
        float_fc = []
        for i in fhd:
            i = i.strip()
            if i.startswith("# Redundant rate in treatment"):
                temp = i.split(":")
                self.redundant_ratio = str(1-float(temp[1]))
            if i and not i.startswith("#") and not i.startswith("chr\t"):
                fs = i.split("\t")
                fc = fs[7]
                float_fc.append(float(fc))
        d = sorted(float_fc)
        d20 = [x for x in d if x >= 20]
        d10 = [x for x in d if x >= 10]
        self.totalpeaks = len(d)
        print d20[1:10]
        self.fold_20 = len(d20)
        self.fold_10 = len(d10)
        print self.fold_20
        print self.fold_10
        print self.totalpeaks
        print self.redundant_ratio
        fhd.close()

    def _high_confidentPeaks_info(self):
        """
        """
        
        historyDataName = os.path.split(chilin.__file__)[0] + '/' + 'db/lg_fold_10.txt'
        pdfName = self.path['qcresult']['fold_ratio']
        rCode = self.path['qcresult']['fold_ratio_r']
        lg_10 = math.log(self.fold_10,10)
        historyData = open(historyDataName).readlines()[0].strip()
        f = open(rCode,'w')
        f.write('peaks_fc <- c(%s)\n' %historyData)
        f.write('fn <- ecdf(peaks_fc)\n')
        f.write('density <- fn(peaks_fc)\n')
        f.write("pdf('%s')\n" %pdfName)
        f.write("plot(ecdf(peaks_fc),verticals=TRUE,col.hor='blue', col.vert='black',main='Fold 10 peaks Distribution',xlab='lg(number of fold_enrichment>10 peaks)',ylab='Cumulative density function of all public data')\n")
        f.write("points(%f,fn(%f),pch=21,bg='red')\n" % (lg_10,lg_10))
        f.write('abline(v=3,col="red")\n')
        f.write("text(3.5,0,'cutoff=3')\n")
        f.write('dev.off()\n')
        f.close()
        call('Rscript %s' % rCode, shell = True)
        

    def _velcro_ratio_info(self,peakbed):
        """verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
         The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""

        historyFile = os.path.split(chilin.__file__)[0] + '/' + 'db/wgEncodeHg19ConsensusSignalArtifactRegions.bed'
        overlapped_bed_file = "overlapped_bed_file"
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
        rCode = self.path['qcresult']['velcro_ratio_r']
        pdfName = self.path['qcresult']['velcro_ratio']
        historyData = self.historyData[3][0:-1]
        f = open(rCode,'w')
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n"% pdfName)
        f.write("rawdata<-c(%s)\n" % historyData)
        f.write("fn<-ecdf(rawdata)\n")
        f.write("plot(ecdf(rawdata), verticals=TRUE,col.hor='blue', col.vert='black',main='velro ratio',xlab='velcro ratio',ylab='Fn(velcro ratio)')\n")
        f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(velcro_ratio,velcro_ratio,int(pch[0]),col[0]))
        f.write("dev.off()\n")
        f.close()
        call('Rscript %s' % rCode, shell = True)
        



    def _DHS_ratio_info(self,peakbed):
        """ DHS ratio indicate the percentage of peaks overlap with DHSs site.
        The function can describe  the particularly dataset's DHS ratio quality of all historic data.
        """
        historyFile = os.path.split(chilin.__file__)[0] + '/' + 'db/DHS_hg19.bed' 
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
        print num_overlapped_peaks
        print dhs_ratio
        historyData = self.historyData[2][0:-1]
        rCode = self.path['qcresult']['dhs_ratio_r']
        pdfName = self.path['qcresult']['dhs_ratio']
        f = open(rCode,'w')
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        f.write("pdf('%s',height=8.5,width=8.5)\n" % pdfName)
        f.write("rawdata<-c(%s)\n" % historyData)
        f.write("fn<-ecdf(rawdata)\n")
        f.write("plot(ecdf(rawdata), verticals=TRUE,col.hor='blue', col.vert='black',main='overlapped_with_DHSs',xlab='overlapped_with_DHSs',ylab='Fn(overlapped_with_DHSs)')"+"\n")
        f.write("points(%f,fn(%f),pch=%d,bg='%s')\n" %(dhs_ratio,dhs_ratio,int(pch[0]),col[0]))
        f.write("dev.off()\n")
        f.close()
        call('Rscript %s' % rCode, shell = True)
    def _replicate_info(self,vennGraph = '',correlationPlot = ''):
        """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""
        self.replicte_check = True
        self.replicate_venn_stat = vennGraph
        self.replicate_corr_stat = correlationPlot
        print 'replicate_info\n'
    def _render(self):
        rend = {}
        rend['PeakcallingQC_check'] = True
        rend['hight_confident_peak_graph'] = 'adasfsdfadsfs'
        temp = self.template.render(rend)

    def run(self,peaksxls,peaksbed):
        """ Run some PeakcallingQC function to get final result. """

        historyDataName = os.path.split(chilin.__file__)[0] + '/' + 'db/all_data.txt'
        self.historyData = open(historyDataName).readlines()
        self._peak_summary_info(peaksxls)
        self._high_confidentPeaks_info()
        self._velcro_ratio_info(peaksbed)
        self._DHS_ratio_info(peaksbed)
        self._replicate_info()
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
        self.texfile = texfile
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
        
        # plot 
        rCode = 'temp.r'
        pdfName = 'ceasa.pdf'
        f = open(rCode,'w')
        list_fcr = ','.join(list_fc)
        f.write("pdf('test.pdf',height=11.5,width=8.5)\n")
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

        f.write("pdf('dis.pdf',height=7,width=5.5)\n")
        f.write("nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), width= c(1,1),height=c(1,1),respect=TRUE)\n")
        Metascript = '\n'.join(Metaregxcontent.split('\n')[1:]) + '\n'
        f.write(Metascript)
        f.write('dev.off()\n')
        call("Rscript %s"% rCode, shell = True)
        
        f.close()
        
        print 'ceas qc'
    def _conservation_info(self):
        """ Density plot of peaks's conservation."""
        print 'conservation qc'

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
        outdir = self.path[qcresult][folder]
        species = 'Homo sapiens'
        zipFile = zipfile.ZipFile(Zippath)
        data = zipFile.read('results/mdseqpos_out.html')
        if  not os.path.isdir(outdir+"seqposimg/"):
            cmd = 'unzip ' + Zippath + ' -d ' + outdir+'seqposimg' + ' \"results/img/*\"'
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
                        output.append([i['factors'][0].upper(),str(i['hits']),str(i['zscore']), i['id']])
                        if len(output) >= 7:
                            break
        print output
        print 'motif info\n'
        print 'motif info\n'
    def _render(self):
        pass
    def run(self):
        """ Run some AnnotationQC function. """
        peaksxls = '/Users/Samleo/mybin/testdata/1277_peaks.xls'
        ceasCode = '/Users/Samleo/mybin/testdata/1277_ceas.R'
        Zippath = '/Users/Samleo/mybin/testdata/1277_seqpos.zip'
        self._ceas_info(peaksxls,ceasCode)
        self._conservation_info()
        self._motif_info(Zippath)
    def check(self):
        """ Check whether AnnotationQC's result is ok. """
        print 'pass or not'


