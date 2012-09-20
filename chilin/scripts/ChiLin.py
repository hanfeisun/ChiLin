#!/usr/bin/env python
"""
Main program of DC pipeline
"""
import sys
import os
from pkg_resources import resource_filename
import argparse
from functools import partial
from chilin.dc import (gen_conf,
                       PipeGroom,
                       PipePreparation,
                       PipeBowtie,
                       PipeMACS2,
                       PipeVennCor,
                       PipeCEAS,
                       PipeConserv,
                       PipeMotif,
                       package
		       )

from chilin.qc import (RawQC,
                       MappingQC,
                       PeakcallingQC,
                       AnnotationQC,
                       SummaryQC)

class ChiLinParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()


def parse_args():
    description = "ChiLin : A clear ChIP-seq pipeline"
    parser = ChiLinParser(description = description)
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    template_parser = sub_parsers.add_parser("gen",  help = "generate a template of config file",
                                             description = "ChiLin-gen: A config template generator for ChiLin")
    template_parser.add_argument("-s","--species", choices = ("hg19", "mm9"), required = True)

    pipe_parser = sub_parsers.add_parser("run", help = "run pipeline using a config file",
                                         description = "ChiLin-run: Run ChiLin pipeline using a config file")
    pipe_parser.add_argument("-c","--config", required = True,
                        help = "specify the config file to use", )
    pipe_parser.add_argument("-t", dest = "atype", choices = ("Dnase", "Histone", "TF"), required = True,
                        help = "the most important option for ChiLin specify the analysis type and the shiftsize {Dnase: 50, Histone and TF:73} of MACS2")
    pipe_parser.add_argument("-m", dest = "cor_method", default = "mean",choices = ("mean",),
                        help = "specify method for correlation plot")
    pipe_parser.add_argument("-e", dest = "model", action = "store_true", default = False,
                             help = "MACS2 establish model or not, use default shift-size if failed.")
    pipe_parser.add_argument("-p", dest = "top_peaks", type = int, default = 5000,
                        help = "specify peaks number for CEAS")
    pipe_parser.add_argument("--threads", dest = "max_threads", type = int, default = 1, choices = range(1,9),
                             help = "How many threads can be used")

    pipe_parser.add_argument("-d", dest = "step_end", type = int, default = 6,
                        help = "specify the end step of pipeline, 1 for bowtie, 2 for macs, 3 for venn and correlation, 4 for ceas, 5 for conservation, 6 for motif, Note: if you only have bed file, start from 2")
    pipe_parser.add_argument("--debug", help = "debug mode", action = "store_true", default = False)
    
    args = parser.parse_args()


    if args.sub_command == "gen":
        gen_conf(args.species)
        sys.exit(0)

    if args.sub_command == "run":
        shift = {"TF": 73, "Histone": 73, "Dnase": 50}
        args.shiftsize = shift[args.atype]
    return args

def main():
    args = parse_args()
    
    pp = PipePreparation(args.config)
    log_func = pp.func_log()

    pp.check_conf()

    conf = pp.get_config()
    rule = pp.get_rule()

    
    p = lambda func:partial(func, conf=conf, rule=rule, log=log_func, debug=args.debug,
                            texfile=pp.get_tex(),
                            datasummary = pp.get_summary(),
                            stepcontrol=args.step_end,
                            a_type = args.atype,
                            shiftsize = args.shiftsize,
                            peaksnumber = args.top_peaks,
                            ArgsionMethod = args.cor_method,
                            Macs2Model = args.model,
                            threads = args.max_threads)
    print conf
    groom = p(PipeGroom)()
    groom.run()

    rawqc = p(RawQC)()
    rawqc.run()

    bowtie = p(PipeBowtie)()
    bowtie.run()

    mappingqc = p(MappingQC)(summarycheck = rawqc.summarycheck)
    mappingqc.run()

    macs2 = p(PipeMACS2)()
    macs2.run()

    pipevenncor = p(PipeVennCor)(ratios = macs2.rendercontent)
    pipevenncor.run()

    peakcallingqc = p(PeakcallingQC)(summarycheck = mappingqc.summarycheck)
    peakcallingqc.run()
    
    pipeceas = p(PipeCEAS)()
    pipeceas.run()

    pipeconserv = p(PipeConserv)()
    pipeconserv.run()

    pipemotif = p(PipeMotif)()
    pipemotif.run()

    annotationqc = p(AnnotationQC)(summarycheck = peakcallingqc.summarycheck)
    annotationqc.run(args.atype)


    summaryqc = p(SummaryQC)()
    summaryqc.run(annotationqc.summarycheck)
    p(package)()
            


if __name__ == "__main__":
    try: 
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interupt~Bye")
        sys.exit(0)
