#!/usr/bin/python2.6
# Time-stamp: <2011-03-30 13:58:29 Tao Liu>

"""Description: Draw correlation plot for many wiggle files.

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import re
import logging
from optparse import OptionParser
import urllib2
import tempfile
import gzip
import subprocess
from CistromeAP.taolib.CoreLib.Parser import WiggleIO
from CistromeAP.taolib.CoreLib.BasicStat.Func import * 
import bisect
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )


UCSC_chrom_URL = r'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/chromInfo.txt.gz'

# ------------------------------------
# Misc functions
# ------------------------------------

error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info

def get_chrom_length ( dbname ):
    # first try to find a local file called dbname
    try:
        fhd = open(dbname,"r")
        chrom_len = {}
        for l in fhd:
            fs = l.split()
            try:
                chrom_len[fs[0]] = int(fs[1])
            except:
                pass
        fhd.close()
    except:
        # get chromosome length from UCSC db download page
        f = urllib2.urlopen(UCSC_chrom_URL % (dbname))
        tmpfname = tempfile.mkstemp(prefix="qcmanychip")[1]
        tmpf = open(tmpfname,'w')
        tmpf.write(f.read())                # write file content to temp file
        tmpf.close()
        f.close
        # read it
        fhd = gzip.open(tmpfname,'r')
        chrom_len = {}
        for l in fhd:
            fs = l.split()
            try:
                chrom_len[fs[0]] = int(fs[1])
            except:
                pass
        fhd.close()
        os.unlink(tmpfname)
    return chrom_len

def Binned(step_region, method, index):
    step = index[1]-index[0]
    result = []
    region_starts = [t[0] for t in step_region]
    region_value = [t[2] for t in step_region]
    length = len(region_value)
    for i in index:
        left = bisect.bisect_right(region_starts, i)
        right = bisect.bisect_right(region_starts, i+step)
        if right == length:
            sub = region_value[left:]
        elif left == length:
            sub = []
        else:
            sub = region_value[left:right]
        if sub:
            result.append(str(method(sub)))
        else:
            result.append("NA")
    return result
    
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog <-r rfile> [options] <wiggle files> ..."
    description = "Draw correlation plot for many wiggle files. Based on qc_chIP_whole.py"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-d","--db",type="str",dest="dbname",help="UCSC db name for the assembly. Default: ce4",default="ce4")
    optparser.add_option("-r","--rfile",dest="rfile",
                         help="R output file. If not set, do not save R file.")
    optparser.add_option("-s","--step",dest="step",type="int",
                         help="sampling step in kbps. default: 100, minimal: 1",default=100)
    optparser.add_option("-z","--imgsize",dest="imgsize",type="int",
                         help="image size in inches, note the PNG dpi is 72. default: 10, minimal: 10",default=10)    
    optparser.add_option("-f","--format",dest="imgformat",type="string",
                         help="image format. PDF or PNG",default='PDF')
    optparser.add_option("-m","--method",dest="method",type="string",default="median",
                         help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean', and 'sample' (just take one point out of a data set). Default: median")
    optparser.add_option("-l","--wig-label",dest="wiglabel",type="string",action="append",
                         help="the wiggle file labels in the figure. No space is allowed. This option should be used same times as wiggle files, and please input them in the same order as -w option. default: will use the wiggle file filename as labels.")
    optparser.add_option("--min-score",dest="minscore",type="float",default=-10000,
                         help="minimum score included in calculation. Points w/ score lower than this will be discarded.")
    optparser.add_option("--max-score",dest="maxscore",type="float",default=10000,
                         help="maximum score included in calculation. Points w/ score larger than this will be discarded.")
    optparser.add_option("-H","--heatmap",dest="heatmap",action="store_true",default=False,
                         help="If True, a heatmap image will be generated instead of paired scatterplot image.")
    
    (options,wigfiles) = optparser.parse_args()

    imgfmt = options.imgformat.upper()
    if imgfmt != 'PDF' and imgfmt != 'PNG':
        print "unrecognized format: %s" % imgfmt
        sys.exit(1)

    method = options.method.lower()
    if method == 'median':
        medfunc = median
    elif method == 'mean':
        medfunc = mean
    elif method  == 'sample':
        medfunc = lambda u: u[-1]
    else:
        print "unrecognized method: %s" % (method)
        sys.exit(1)

    wigfilenum = len(wigfiles)
    if wigfilenum < 2 or not options.rfile:
        error("must provide >=2 wiggle files")
        optparser.print_help()
        sys.exit(1)

    # wig labels
    if options.wiglabel and len(options.wiglabel) == wigfilenum:
        wiglabel = options.wiglabel
    else:  # or use the filename
        wiglabel = map(lambda x:os.path.basename(x),wigfiles)
        
    if options.step < 1:
        error("Step can not be lower than 1!")
        sys.exit(1)
    if options.imgsize < 10:
        error("Image size can not be lower than 10!")
        sys.exit(1)

    # check the files
    for f in wigfiles:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)
        
    #wigfhds = map(open,wigfiles)        # file handlers for wiggle files
    info("number of wiggle files: %d" % wigfilenum)
    # get chromosome length info from UCSC
    info("connect to UCSC to get chromosome length information")
    try:
        chrom_len = get_chrom_length(options.dbname)
    except:
        error("Error!")
        sys.exit(1)

    chromsdict = {}
    for iwig in wigfiles:
        p = subprocess.Popen(['bigWigInfo','-chroms', iwig],stdout=subprocess.PIPE).communicate()[0].split('\n')
        ichroms = [t.split()[0] for t in p if len(t.split())==3 and t.split()[0].startswith('chr')]
        for i in ichroms:
            chromsdict[i] = chromsdict.setdefault(i,0)+1
    
    chroms = [] #store common chroms in every wig file.
    for c in chromsdict.keys():
        if chromsdict[c]==wigfilenum:
            chroms.append(c)
    if not chroms:
        error('No common chrom found')
        sys.exit()
    info("common chromosomes are %s." % ",".join(chroms))

    if options.rfile:
        rfhd = open(options.rfile,"w")
        rfhd.write('''require("RColorBrewer") ## from CRAN\n''')

    # for each wig file, sample...
    for i in range(len(wigfiles)):
        #wigfhd = wigfhds[i]
        #wigfhd.seek(0)                  # reset
        bw = BigWigFile(open(wigfiles[i]))
        
        info("read wiggle track from bigwig file #%d" % (i+1))
        #wig = WiggleIO.WiggleIO(wigfhd).build_binKeeper(chromLenDict=chrom_len,binsize=options.step*1000)
        profiles = []  # size of p is number of cages in binkeeper
        for chrom in chroms:
            try:
                end_of_chrom = chrom_len[chrom] +10000000
            except:
                warn("chromosome %s can not be found in UCSC... skip..." % chrom)
                continue
            step_region = bw.get(chrom, 0, end_of_chrom)
            index = range(0, end_of_chrom - options.step*1000, options.step*1000)
            profiles.append(Binned(step_region, medfunc, index)) #list include "NA" and str num.
            
        info("write values to r file")
        rfhd.write("p%d <- c(%s)\n" %(i, ','.join(profiles[-1])))
        
    rfhd.write("c <- cbind(p0")
    for i in range(wigfilenum-1):
        rfhd.write(",p%d" % (i+1))
    rfhd.write(")\n")
    
    rfhd.write("c <- c[ c[,1]<=%f & c[,1]>=%f " % (options.maxscore,options.minscore))
    for i in range(wigfilenum-1):
        rfhd.write("& c[,%d]<=%f & c[,%d]>=%f " % (i+2,options.maxscore,i+2,options.minscore))
    rfhd.write(",]\n")
    if imgfmt == 'PDF':
        rfhd.write("pdf(\"%s.pdf\",width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))
    elif imgfmt == 'PNG':
        rfhd.write("png(\"%s.png\",units=\"in\",res=150,width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))

    if options.heatmap:                 # heatmap
        rfhd.write('library(gplots)\n')
        rfhd.write('''
m <- cor(c, method="pearson", use="pairwise.complete.obs")
''')
        labels = ",".join(map(lambda x:"\""+x+"\"",wiglabel))
        rfhd.write("rownames(m) <- c(%s)\n" % labels)
        rfhd.write("colnames(m) <- c(%s)\n" % labels)         
        rfhd.write('# draw the heatmap using gplots heatmap.2\n') 
#        rfhd.write('bitmap("%s.bmp",width=%d,height=%d)\n' % (options.rfile,options.imgsize,options.imgsize))
        rfhd.write('mn <- -1\n')
        rfhd.write('mx <- 1\n')
        rfhd.write('n <- 98\n')
        rfhd.write('bias <- 1\n')
        rfhd.write('mc <- matrix(as.character(round(m, 2)), ncol=dim(m)[2])\n')
        rfhd.write('breaks <- seq(mn, mx, (mx-mn)/(n))\n')
        rfhd.write('cr <- colorRampPalette(colors = c("#2927FF","#FFFFFF","#DF5C5C"), bias=bias)\n')
        rfhd.write('heatmap.2(m, col = cr(n), breaks=breaks, trace="none", cellnote=mc, notecol="black", notecex=1.8, keysize=0.5, density.info="histogram", margins=c(27.0,27.0), cexRow=2.20, cexCol=2.20, revC=T, symm=T)\n')
    else:                               # scatterplot
        rfhd.write('''
panel.plot <- function( x,y, ... )
{
  par(new=TRUE)
  m <- cbind(x,y)
  plot(m,col=densCols(m),pch=20)
  lines(lowess(m[!is.na(m[,1])&!is.na(m[,2]),]),col="red")  
}
    
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,use="complete.obs")
  txt <- format(round(r,2),width=5,nsmall=2)
  #format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  text(0.5, 0.5, txt, cex = cex.cor)
}
''')
        labels = ",".join(map(lambda x:"\""+x+"\"",wiglabel))
        rfhd.write('''
pairs(c, lower.panel=panel.plot, upper.panel=panel.cor, labels=c(%s))
''' % (labels))

    rfhd.write("dev.off()\n")
    rfhd.close()

    # try to call R
    try:
        subprocess.call(['Rscript',options.rfile])
    except:
        info("Please check %s" % options.rfile)
    else:
        info("Please check %s" % (options.rfile+'.bmp'))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
