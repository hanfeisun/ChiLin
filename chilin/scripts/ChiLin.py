#!/usr/bin/python
"""
Main program of DC pipeline
"""
from chilin.dc import *
from chilin.qc import *
from optparse import OptionParser
from subprocess import call

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
    Preparation.checkconf()
    conf = Preparation.ChiLinconfi
    PathFinder(conf['userinfo']['treatpath'],conf['userinfo']['controlpath'])
    outputd = conf['userinfo']['outputdirectory']
    print outputd
    call('mkdir %s & cd %s' % (outputd, outputd), shell = True)
    texfile = open('tex.tex', 'wb')

#    log = LogWriter(open('log', 'w'))
#    judge = Preparation.checkconf()
##    if judge == False:
##        sys.exit()
#    fastqc_check = RawQC(texfile).run()
    print Preparation.ChiLinconfigs
    RawQC(texfile).run()
    texfile.close()



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt():
        print "User stops me:)"
    finally:
        print "Welcome to ChiLin"

@
