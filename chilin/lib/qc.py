import chilin
from chilin.dc import *
import os
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
    def _infile_parse(self,dataname):
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
            fastqc_out = self.path['qcresult']['folder']+'/'+d.split('/')[-1]+'_fastqc'
            changed_name = self.path['qcresult']['folder']+'/'+names[i]+'_fastqc'
            cmd = 'mv {0} {1}'
            cmd = cmd.format(fastqc_out,changed_name)
            call(cmd,shell=True)
            call('rm %s.zip'% fastqc_out,shell=True)
            dataname = changed_name+'/fastqc_data.txt'
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
        historyData = os.path.split(chilin.__file__)[0] + '/' + 'db/fastqc_value_list.txt'
        inf=open(historyData,'rU')
        peaklist1=inf.readline()
        oufe=open(rCode,'w')
        oufe.write("setwd('%s')\n" %self.conf['userinfo']['outputdirectory'])
        oufe.write("sequence_quality_score<-c(%s)\n" % str(peaklist1)[1:-1])
        col=['#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054']
        pch=[21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25]
        oufe.write("ecdf(sequence_quality_score)->fn\n")
        oufe.write("fn(sequence_quality_score)->density\n")
        oufe.write("cbind(sequence_quality_score,density)->fndd\n")
        oufe.write("fndd2<-fndd[order(fndd[,1]),]\n")
        oufe.write("pdf('%s')\n" % pdfName)
        oufe.write("plot(fndd2,type='b',pch=18,col=2,main='Sequence Quality Score Cumulative Percentage',ylab='cummulative density function of all public data')\n")
        j=0
        for p in npeakl:
            oufe.write("points(%d,fn(%d),pch=%d,bg='%s')\n" %(int(p),int(p),int(pch[j]),col[j]))
            j=j+1
        oufe.write("legend('topleft',c(%s),pch=c(%s),pt.bg=c(%s))\n" %(str(names)[1:-1],str(pch[:len(names)])[1:-1],str(col[:len(names)])[1:-1]))
        oufe.write("dev.off()\n")
        oufe.close()
        inf.close()
        call('Rscript %s' % pdfName, shell = True)
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

    def _basic_mapping_statistics_info(self):

        """ Stastic summary of mapping result for each sample. """
        print 'basic_mapping_statistics'
        self.mappable_summary_stat = 'basic_mapping_table'
        return self.mappable_summary_stat

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
                mappable_ratio_graph = self.mappable_ratio_stat, redundant_ratio_graph = self.redundant_ratio_stat)
        self.filehandle.write(temp)
        self.filehandle.flush()
    def run(self):
        
        self.MappingQC_check = True
        names = ['rep1','rep2','rep3']
        mapped_ratio = [0.9,0.8,0.8]
        redundant_ratio = [0.1,0.2,0.2]
        historyData = os.path.split(chilin.__file__)[0] + '/' + 'db/all_data.txt'
        self.historyData = open(historyData).readlines()
        """ Run some MappingQC function to get final result. """
        self.mappable_summary_stat = self._basic_mapping_statistics_info()
        self.mappable_ratio_stat = self._mappable_ratio_info(mapped_ratio,names)
        self.redundant_ratio_stat = self._redundant_ratio_info(redundant_ratio,names)
        self._render()
    def check():
        """Check whether the MappingQC's result is ok. """
        print 'mapping qc pass or not'


class PeakcallingQC(QC_Controller):
    """ PeakcallingQC aims to describe the quality of peak calling result."""
    def __init__(self,configs = ''):
        super(PeakcallingQC, self).__init__()
        self.conf = configs
        print 'init peak calling  qc'
    def _peak_summary_info(self):
        """Basic statistic of peak calling result."""
        print ' summary_info'

    def _velcro_ratio_info(self):
        """verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
         The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""
        print 'velcro_ratio_info '

    def _DHS_ratio_info(self):
        """ DHS ratio indicate the percentage of peaks overlap with DHSs site.
        The function can describe  the particularly dataset's DHS ratio quality of all historic data.
        """
        print 'DHS_ratio_info'
    def _replicate_info(self):
        """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""
        print 'replicate_info\n'
    def _render(self):
        pass

    def run(self):
        """ Run some PeakcallingQC function to get final result. """
        self._peak_summary_info()
        self._velcro_ratio_info()
        self._DHS_ratio_info()
        self._replicate_info()
    def check():
        """ Check whether PeakcallingQC's result is ok. """
        print 'pass or not'

        
class AnnotationQC(QC_Controller):
    """ AnnotationQC aims to describe the quality of annotations after peak calling. """ 
    def __init__(self,configs = ''):
        super(AnnotationQC, self).__init__()
        self.conf = configs
        print 'intialization of function qc'
    def _ceas_info(self):
        """ Describe peaks' distribution and relative position. """
        print 'ceas qc'
    def _conservation_info(self):
        """ Density plot of peaks's conservation."""
        print 'conservation qc'
    def _motif_info(self):
        """ QC of Sepose. """
        print 'motif info\n'
    def _render(self):
        pass
    def run(self):
        """ Run some AnnotationQC function. """
        self._ceas_info()
        self._conservation_info()
        self._motif_info()
    def check(self):
        """ Check whether AnnotationQC's result is ok. """
        print 'pass or not'


