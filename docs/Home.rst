=================
Home
=================

Project Purpose
===================

* ChIP-seq experiment has been a mature and wide-spread technique for detecting the TF and histone modification distribution from the genome scale.
* Along with the popularity of the technique and the increasingly huge number of highthroughput datasets, it may be confusing for biologists to get a quick and easy access to understand their biological meaning, and the same important thing is the unbiased judgement of the data quality.So, it's necessary for us establishing a universal and easy-to-use ChIP-seq data analysis pipeline for biologists.
* For bioinformaticians, you may find this may be the most extensible and flexible ChIP-seq integrated python packages ever. It provides various genomic data format support and high-rated analysis toolbox. Any bug report is welcome.


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
.. _ChiLin:: https://bitbucket.org/shenglinmei/chilin
