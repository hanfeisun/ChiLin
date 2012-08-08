#!/usr/bin/env python
"""
Main program of DC pipeline
"""
import sys
import os
from pkg_resources import resource_filename
from optparse import OptionParser

from chilin.dc import (PipePreparation,
                       PipeBowtie,
                       PipeMACS2,
                       PipeVennCor,
                       PipeCEAS,
                       PipeConserv,
                       PipeMotif,
                       package)

from chilin.qc import (RawQC,
                       MappingQC,
                       PeakcallingQC,
                       AnnotationQC,
                       SummaryQC)

def parse_options():
    usage = "usage: %prog <ChiLin.conf Path> [optional]"
    description = "ChiLin : A clear ChIP-seq pipeline"
    parser = OptionParser(version = "%prog 1.0.0", description = description, usage = usage)
    parser.add_option("-m", dest = "cor_method", type = "string", default = "mean",
                      help = "specify method for correlation plot")
    parser.add_option("-p", dest = "top_peaks", type = "string", default = '5000',
                      help = "specify peaks number for CEAS")
    parser.add_option("-t", dest = "atype", type = "string",
                      help = "It is a most important option for ChiLin specify the analysis type, supported Dnase, Histone, TF")
    parser.add_option("-s", dest = "shiftsize", type = "string", default= '73',
                      help = "specify the fixed shiftsize for MACS2, advice Dnase: 50, Histone and TF:73")
    parser.add_option("-d", dest = "step_end", type = "int", default = 6,
                      help = "specify the end step of pipeline, 1 for bowtie, 2 for macs, 3 for venn and correlation, 4 for ceas, 5 for conservation, 6 for motif, Note: if you only have bed file, start from 2")
    parser.add_option("-e", action="store_true", dest="example_conf", default=False,
                      help = "Show an example of configuration file")
    (options, args) = parser.parse_args()

    if options.example_conf:
        with open(resource_filename("chilin", os.path.join("db", "ChiLin_human.conf"))) as f:
            print f.read()
        sys.exit(1)
    if not args:
        parser.print_help()
        sys.exit(1)
    if not options.atype :
        parser.print_help()
        print "Please select the analysis type (Dnase, Histone or TF)!"
        sys.exit(1)
        
    if options.atype in ['TF', 'Histone']:
        options.shiftsize = '73'
    elif options.atype == 'Dnase':
        options.shiftsize = '50'
    return options, args

def main():
    opt, args = parse_options()
    conf_file = args[0]
    
    pp = PipePreparation(conf_file)
    log_func = pp.func_log()

    pp.check_conf()
    conf = pp.get_config()
    rule = pp.get_rule()

    with open(pp.get_tex(),"w") as f_tex:
        with open(pp.get_summary(),"w") as f_sum:
            rawqc = RawQC(conf = conf, rule = rule,
                          texfile = f_tex, log = log_func)
            rawqc.run()

            pipebowtie = PipeBowtie(conf = conf, rule = rule, log = log_func,
                                datasummary= f_sum, stepcontrol = opt.step_end)
            pipebowtie.run()

            mappingqc = MappingQC(conf = conf, rule = rule, log = log_func,
                                  texfile = f_tex, summarycheck = rawqc.summarycheck)
            mappingqc.run(bedft)
            
            macs2 = PipeMACS2(conf = conf, rule = rule, log = log_func,
                              datasummary = f_sum, stepcontrol = opt.step_end,
                              shiftsize = opt.shiftsize, bedft = bedft)
            macs2.run()
            
            pipevenncor = PipeVennCor(conf = conf, rule = rule, log = log_func,
                                      datasummary = f_sum, stepcontrol = opt.step_end,
                                      ratios = macs2.rendercontent, peaksnumber = opt.top_peaks,
                                      OptionMethod = opt.cor_method)
            pipevenncor.run()
            
            peakcallingqc = PeakcallingQC(conf = conf,rule =rule, log = log_func,
                                          texfile = f_tex, summarycheck = mappingqc.sum_check)
            peakcallingqc.run()
            
            pipeceas = PipeCEAS(conf = conf, rule = rule, log = log_func,
                                stepcontrol = opt.step_end, peaksnumber = p)
            pipeceas.run()
            
            pipeconserv = PipeConserv(conf = conf, rule = rule, log = log_func,
                                      stepcontrol = opt.step_end, a_type = opt.atype)
            pipeconserv.run()
            
            pipemotif = PipeMotif(conf = conf, rule = rule, log = log_func,
                                  stepcontrol = opt.step_end)
            pipemotif.run()
            
            annotationqc = AnnotationQC(conf = conf, rule = rule, log = log_func,
                                        texfile = f_tex, summarycheck = peakcallingqc.sum_check)
            annotationqc.run()
            
            summaryqc = SummaryQC(conf = conf, rule = rule, log = log_func,
                                  texfile = f_tex)
            summaryqc.run(sum_check)
            


if __name__ == "__main__":
    try: 
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interupt~Bye")
        sys.exit(0)


