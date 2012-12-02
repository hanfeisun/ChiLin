=================
Home
=================

Project Purpose
===================

* ChIP-seq experiment has been a mature and wide-spread technique for detecting the TF, histone modification and chromatin factors distribution from the genome scale.
* Along with the popularity of the technique and the increasingly huge number of highthroughput datasets, it may be confusing for biologists to get a quick and easy access to understand their biological meaning, and the same important thing is the unbiased judgement of the data quality.So, it's necessary for us to establish a universal and user-friendly ChIP-seq data analysis pipeline for biologists.
* For bioinformaticians, this is the most extensible and flexible ChIP-seq implemented with python so far. It support various genomic data format and integrate high-rated analysis toolbox.

Project Goal
==================
There are two main command for the pipeline.
The one `merge_DC` is for previous reservation version for constructing Cistrome Secondary Database;
The New Version  ChiLin_ provide a more flexible handle for understanding the ChIP-seq analysis workflow.

There are two layers for pipeline goal.This is part of the Cistrome_ project

1.command line ::

  The goal is to simply input the sequence data files up to input format
  support and fill in the customized table of experiment descriptions,
  output the desired format of DC and QC report

2.web server ::

  The ChIP-seq pipeline could be incorporated into the Cistrome

.. _Cistrome: http://Cistrome.org
.. _ChiLin: https://bitbucket.org/shenglinmei/chilin
