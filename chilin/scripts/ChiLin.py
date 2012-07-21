#!/usr/bin/python
"""
Main program of DC pipeline
"""
from chilin.dc import *
from chilin.qc import *
from optparse import OptionParser
from subprocess import call
import os

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
                      help = "specify the fixed shiftsize for MACS2")
#    parser.add_option("-s", dest = "stepcontrol", type = "store_true", default = True)
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
    outputd = conf['userinfo']['outputdirectory']
    print outputd, 'here'
    if not os.path.exists(outputd):
        call('mkdir %s & cd %s' % (outputd, outputd), shell = True)
    else:
        call('cd %s' % outputd, shell = True)

#   log = LogWriter('log')
#    log.record('test')
    Path = PathFinder(outputd, conf['userinfo']['datasetid'], conf['userinfo']['treatpath'], conf['userinfo']['controlpath'])
    fastqcname = Path.qcfilepath()
    print fastqcname

    texfile = open('tex.tex', 'wb')


    judge = Preparation.checkconf()
    print "Test for output"
    if judge == False:
        sys.exit()
    print outputd
    fastqc_check = RawQC(texfile).run(conf['qc'],fastqcname,conf['userinfo']['treatpath'], conf['userinfo']['controlpath'],conf['userinfo']['outputdirectory'])
    print Preparation.ChiLinconfigs
    texfile.close()

    fastqc_judge = True
    if fastqc_judge == True:
        bowtiename = Path.bowtiefilepath()
        print bowtiename
        PipeBowtie().summary(
        

    



if __name__ == '__main__':
    print "Welcome to ChiLin"
    try:
        main()
    except KeyboardInterrupt():
        print "User stops me:)"
