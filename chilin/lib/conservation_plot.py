#!/usr/bin/python
# Time-stamp: <2011-07-15 21:12:20 Jian Ma>

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
import subprocess
import math
from optparse import OptionParser

from CistromeAP.taolib.CoreLib.Parser import WiggleIO, BedIO
from CistromeAP.taolib.CoreLib.BasicStat.Func import * 

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
#bigWigSummary = 'bigWigSummary'

# ------------------------------------
# Misc functions
# ------------------------------------

error  = logging.critical		# function alias
warn   = logging.warning
debug  = logging.debug
info   = logging.info

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog <-d path> [options] <bed files> ..."
    description = "Draw conservation plot for many bed files."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option('-H','--height', dest='height',type='int',default=10, help="height of plot")
    optparser.add_option('-W','--width',dest='width',type='int',default=10, help="width of plot")
    optparser.add_option('-w',dest='w',type='int',default=1000, help="window width centered at middle of bed regions,default: 1000")    
    optparser.add_option('-t','--title',dest='title',help="title of the figure. Default: 'Average Phastcons around the Center of Sites'",default= 'Average Phastcons around the Center of Sites')
    optparser.add_option('-d','--phasdb',dest='phasdb',help= 'The directory to store phastcons scores in the server')
    optparser.add_option("-l","--bed-label",dest="bedlabel",type="string",action="append",
                         help="the BED file labels in the figure. No space is allowed. This option should be used same times as -w option, and please input them in the same order as BED files. default: will use the BED file filename as labels.")
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    
    (options,bedfiles) = optparser.parse_args()
    options.pf_res = options.w / 100 # get 100 points to plot
    options.w = options.pf_res * 100 # trim

    bedfiles = map(os.path.abspath,bedfiles)
    bedfilenames = map(os.path.basename,bedfiles)

    bedfilenum = len(bedfiles)

    if bedfilenum < 1 or not options.phasdb:
        optparser.print_help()
        sys.exit(1)

    if options.bedlabel and len(options.bedlabel) == bedfilenum:
        bedlabel = options.bedlabel
    else:                               # or use the filename
        bedlabel = map(lambda x:os.path.basename(x),bedfiles)

    if options.height < 10:
        error("Height can not be lower than 10!")
        sys.exit(1)
    if options.width < 10:
        error("Width can not be smaller than 10!")
        sys.exit(1)

    # check the files
    for f in bedfiles:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)

    # check phastcons db
    if not os.path.isdir(options.phasdb):
        error("%s is not valid!" % options.phasdb)
        sys.exit(1)

    # change wd to phastcons db path
    olddir = os.path.abspath('.')
    os.chdir(options.phasdb)

    phas_chrnames = []
    
    files_phasdb = os.listdir('.')
    for file_phasdb in files_phasdb:
        if file_phasdb.endswith('.bw'):
            name = file_phasdb.rstrip('.bw')
            phas_chrnames.append(name)

    if not phas_chrnames:
        error("%s has no valid phastcons db bw files!" % options.phasdb)
        sys.exit(1)
        
    info("number of bed files: %d" % bedfilenum)

    avgValues = []

    # for each bed file
    for f in bedfiles:
        info("extract phastcons scores using %s" % f)
        scores = extract_phastcons(f,phas_chrnames, options.w, options.pf_res)
        avgValues.append(scores)
    makeBmpFile(avgValues,olddir,options.height,options.width,options.w,options.pf_res,options.title,bedlabel)

def extract_phastcons ( bedfile, phas_chrnames, width, pf_res ):
    """Extract phastcons scores from a bed file.

    Return the average scores
    """
    info("read bed file...")
    bfhd = open(bedfile)
    bed = BedIO.parse_BED(bfhd)

    # calculate the middle point of bed regions then extend left and right by 1/2 width
    bchrs = bed.peaks.keys()
    bchrs.sort()

    chrs = []
    for c in phas_chrnames:
        if c in bchrs:
            chrs.append(c)

    sumscores = []
    for chrom in chrs:
        info("processing chromosome: %s" %chrom)
        pchrom = bed.peaks[chrom]
        bw = BigWigFile(open(chrom+'.bw', 'rb'))
        for i in range(len(pchrom)):
            mid = int((pchrom[i][0]+pchrom[i][1])/2)
            left = int(mid - width/2)
            right = int(mid + width/2)
            
            if left < 0:
                left = 0
                right = width
            
            summarize = bw.summarize(chrom, left, right, width/pf_res)
            if not summarize:
                continue
            dat = summarize.sum_data / summarize.valid_count
            #dat = dat.strip().split('\t')
            sumscores.append(dat)
            
    sumscores = map(list, zip(*sumscores))
    sumscores = [[t2 for t2 in t if not math.isnan(t2)] for t in sumscores]
    conscores = [sum(t)/len(t) for t in sumscores]

    return  conscores
        
def makeBmpFile(avgValues, wd, h,w, width, pf_res, title, bedlabel):
    
    #creating R file in which to write the rscript which defines the correlation plot
    #create and save the file in the current working directory

    fileName = os.path.join(wd, 'tmp')
    rFile = open(fileName+'.R','w')
    bmpname = fileName+'.pdf'
    rscript = 'sink(file=file("/dev/null", "w"), type="message")\n'
    rscript += 'sink(file=file("/dev/null", "w"), type="output")\n'    
    rscript += 'pdf("%s",height=%d,width=%d)\n' %(bmpname,h,w)
    xInfo = range(int(-width/2),int(width/2), pf_res)
    rscript += 'x<-c('+','.join(map(str,xInfo[:-1]))+')\n' # throw the last point which may be buggy
    for i in range(len(avgValues)):
        avgscores = avgValues[i]
        tmpname = 'y'+str(i)
        rscript += tmpname+'<-c('+','.join(map(str,avgscores[:-1]))+')\n' # throw the last point which may be buggy

    tmplist = []
    for i in range(len(avgValues)):
        tmplist.append( "y%d" % i )
    
    rscript += "ymax <- max("+ ",".join(tmplist) +")\n"
    rscript += "ymin <- min("+ ",".join(tmplist) +")\n"    
    rscript += "yquart <- (ymax-ymin)/4\n"

    rscript += 'plot(x,y0,type="l",col=rainbow(%d)[1],main=\"%s\",xlab="Distance from the Center (bp)",ylab="Average Phastcons",ylim=c(ymin-yquart,ymax+yquart))\n' % (len(avgValues),title)
    for i in range(1,len(avgValues)):
        rscript += 'lines(x,y'+str(i)+',col=rainbow(%d)[%d])\n' % (len(avgValues),i+1)
    rscript += 'abline(v=0)\n'
    legend_list = map(lambda x:"'"+x+"'", bedlabel)
    rscript += 'legend("topright",c(%s),col=rainbow(%d),lty=c(%s))\n' % (','.join(legend_list),len(avgValues),','.join(['1']*len(avgValues)))
        
    rscript += 'dev.off()\n'
    rFile.write(rscript)
    rFile.close()
    #executing the R file and forming the pdf file
    data = subprocess.call(['Rscript',fileName+'.R'])

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
