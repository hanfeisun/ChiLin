#!/usr/bin/python
"""
Main program of DC pipeline
"""
from chilin.dc import *
from chilin.qc import *
from optparse import OptionParser
from subprocess import call
import os
import sys
def main():
    usage = "usage: %prog <ChiLin.conf Path> [optional]-m cormethod -p peaksnumber"
    description = "ChiLin : A clear ChIP-seq pipeline"
    parser = OptionParser(version = "%prog 1.0.0", description = description, usage = usage)
    parser.add_option("-m", dest = "cormethod", type = "string", default = "mean",
                      help = "specify method for correlation plot")
    parser.add_option("-p", dest = "peaksnumber", type = "int", default = 3000,
                      help = "specify peaks number for CEAS, Conservation and options")
    parser.add_option("-t", dest = "type", type = "string",
                      help = "specify the analysis type, supported Dnase, Histone, TF")
    parser.add_option("-s", dest = "shiftsize", type = "string", default= '73',
                      help = "specify the fixed shiftsize for MACS2, advice Dnase: 50, Histone and TF:73")
    parser.add_option("-c", dest = "stepcontrol", type = "string", default = 'motif')
    (options, args) = parser.parse_args()

    if not args or not options.type:
        parser.print_help()
        sys.exit('options missing')

    ChiLinConf = args[0]

    Preparation = PipePreparation(ChiLinConf)
    checkresult = Preparation.checkconf()
    if not checkresult:
        sys.exit()
    conf = Preparation.ChiLinconfigs
    conf['trepn'] = len(conf['userinfo']['treatpath'].split(','))
    if conf['userinfo']['controlpath'] != '':
        conf['crepn'] = len(conf['userinfo']['controlpath'].split(','))
    else:
        conf['crepn'] = 0
    print conf['trepn'], conf['crepn']
    outputd = conf['userinfo']['outputdirectory']
    if not os.path.exists(outputd):
        call('mkdir %s' % outputd, shell = True)
        os.chdir(outputd)
    else:
        os.chdir(outputd)

#    log = LogWriter('log')
 #   log.record('test')
    Path = PathFinder(outputd, conf['userinfo']['datasetid'], conf['userinfo']['treatpath'], conf['userinfo']['controlpath'])
    Path.parseconfrep()
    paths = Path.Nameconfigs
    if not conf['userinfo']['outputdirectory'].endswith('/'):
        paths['qcresult']['folder'] = conf['userinfo']['outputdirectory']+'/'+paths['qcresult']['folder']
    else:
        paths['qcresult']['folder'] = conf['userinfo']['outputdirectory']+paths['qcresult']['folder']
    if not os.path.exists(paths['qcresult']['folder']):
        call('mkdir %s' % paths['qcresult']['folder'],shell = True)
    print paths['qcresult']['folder'] 
    texfile = open('tex.tex', 'wb')



    judge = Preparation.checkconf()
    if judge == False:
        sys.exit()
    fastqc_check = RawQC(conf,paths,texfile).run()

    fastqc_judge = True
    if fastqc_judge:

        bowtie = PipeBowtie(conf, paths)
        bowtie.process()
        MappingQC(conf,paths,texfile).run(bowtie.bowtieinfo)
        if bowtie.has_run:
           macs = PipeMACS2(conf, paths)
           macs.process(options.shiftsize)
    #       PeakcallingQC(conf,paths,texfile).run(paths['macstmp']['peaks_xls'],paths['macsresult']['reptreat_peaks'])
           if macs.has_run:
               VennCor = PipeVennCor(conf, paths, peaksnumber, cormethod)
               if VennCor.has_run:
                   CEAS = PipeCEAS(conf, paths)
                   if CEAS.has_run:
                       Motif = PipeMotif(conf, paths)
                       if Motif.has_run:
                           pass
    print bowtie.bowtieinfo
    texfile.close()


if __name__ == '__main__':
    print "Welcome to ChiLin"
    try:
        main()
    except KeyboardInterrupt():
        print "User stops me:)"
