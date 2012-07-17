#!/usr/bin/python
"""
Main program of DC pipeline
"""
from chilin.dc import *
from chilin.qc import *
from optparse import OptionParser

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
    (options, args) = parser.parse_args()

    if not args or not options.type:
        parser.print_help()
        sys.exit(1)
    ChiLinConf = args[0]
    Config = Check().ReadConf(ChiLinConf)

    print Config
#    ['bowtie', 'samtools', 'macs', 'bedtools', 'bed2bam', 'ceas', 'conservation', 'correlation', 'venn', 'seqpos', 'QC']

    print "test OK"





if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt():
        print "User stops me:)"
    finally:
        print "Welcome to ChiLin"




