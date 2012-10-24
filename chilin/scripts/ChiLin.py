#!/usr/bin/env python
"""
Main program of DC pipeline
"""
import sys
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
    #    pipe_parser.add_argument("-e", dest = "model", action = "store_true", default = False,
