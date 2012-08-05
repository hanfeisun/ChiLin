#!/usr/bin/python
"""
Main program of DC pipeline
"""
from chilin.dc import *
from chilin.qc import *
from optparse import OptionParser
import sys

def parse_options():
    usage = "usage: %prog <ChiLin.conf Path> [optional]"
    description = "ChiLin : A clear ChIP-seq pipeline"
    parser = OptionParser(version = "%prog 1.0.0", description = description, usage = usage)
    parser.add_option("-m", dest = "cormethod", type = "string", default = "mean",
                      help = "specify method for correlation plot")
    parser.add_option("-p", dest = "peaksnumber", type = "string", default = '5000',
                      help = "specify peaks number for CEAS")
    parser.add_option("-t", dest = "atype", type = "string",
                      help = "It is a most important option for ChiLin specify the analysis type, supported Dnase, Histone, TF")
    parser.add_option("-s", dest = "shiftsize", type = "string", default= '73',
                      help = "specify the fixed shiftsize for MACS2, advice Dnase: 50, Histone and TF:73")
    parser.add_option("-c", dest = "stepcontrol", type = "int", default = 6,
                      help = "specify the end step of pipeline, 1 for bowtie, 2 for macs, 3 for venn and correlation, 4 for ceas, 5 for conservation, 6 for motif, Note: if you only have bed file, start from 2")
    (options, args) = parser.parse_args()
    if not args or not options.atype or options.atype not in ['TF', 'Dnase', 'Histone']:
        parser.print_help()
        sys.exit(1)
    if options.atype in ['TF', 'Histone']:
        options.shiftsize = '73'
    elif options.atype == 'Dnase':
        options.shiftsize = '50'
    return options, args

def main():
    options, args = parse_options()
    ChiLinConf = args[0]
    Preparation = PipePreparation(ChiLinConf)
    Preparation.checkconf()
    conf = Preparation.ChiLinconfigs
    names = Preparation.Nameconfigs
    log = Preparation.log
    s = options.stepcontrol
    p = options.peaksnumber
    m = options.cormethod
    datasummary = Preparation.DataSummary

    texfile,summarycheck = QC_Controller().QCpreparation(names)
    rawqc = RawQC(conf,names,texfile,summarycheck,log)
    rawqc.run()

    PipeBowtie(conf, names, log, datasummary, s).process()
    mappingqc = MappingQC(conf,names,texfile,rawqc.summarycheck,log)
    mappingqc.run()

    macs2 = PipeMACS2(conf, names, log, datasummary, s, options.shiftsize)
    macs2.process()

    PipeVennCor(conf, names, log, datasummary, s, macs2.rendercontent, p, m).process()
    peakcallingqc = PeakcallingQC(conf,names,texfile,mappingqc.summarycheck,log)
    peakcallingqc.run()

    PipeCEAS(conf, names, log, s, p).process()
    PipeConserv(conf, names, options.atype, log, s).process()
    PipeMotif(conf, names, log, s).process()
    print 'ceas'
    annotationqc = AnnotationQC(conf,names,texfile,peakcallingqc.summarycheck,log)
    annotationqc.run()
    SummaryQC(conf,names,texfile).run(summarycheck)

    package(conf, names, log)



if __name__ == "__main__":
    try: 
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interupt~Bye")
        sys.exit(0)


