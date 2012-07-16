#!/usr/bin/python
"""
Main program of DC pipeline
"""
import sys
sys.path.append("../lib/")
from dc import *
from qc import *
from optparse import OptionParser

def main():
    usage = "usage: %prog <Meta.xls Path> [optional]-m cormethod -p peaksnumber"
    description = "ChiLin : A clear ChIP-seq pipeline"
    parser = OptionParser(version = "%prog 1.0.0", description = description, usage = usage)
    parser.add_option("-m", dest = "cormethod", type = "string", default = "mean",
                      help = "specify method for correlation plot")
    parser.add_option("-p", dest = "peaksnumber", type = "int", default = 3000,
                      help = "specify peaks number for CEAS, Conservation and options")
    (options, args) = parser.parse_args()

    if not args:
        parser.print_help()
        sys.exit(1)
    Meta = args[0]
    Config = Check().ReadConf()
    print Config
#    ['bowtie', 'samtools', 'macs', 'bedtools', 'bed2bam', 'ceas', 'conservation', 'correlation', 'venn', 'seqpos', 'QC']

    
    Meta = Check().MetaParse(Meta)
    

    print "test OK"





if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt():
        print "User stops me:)"
    finally:
        print "Welcome to ChiLin"




