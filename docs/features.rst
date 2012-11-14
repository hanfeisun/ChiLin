========
Features
========
This page is intended for:

**Developers**: in order to be sure they are developing the right project that fulfills requirements provided in this document.

**users**: in order to get familiar with the idea of the project and suggest other features that would make it even more functional.


Cpipe should include the following features.

.. _Raw Data:

Supported formats 
=================

Cpipe should support the following format of ``Raw Data`` as input:

======  ======  ======================
Format  Type    Steps
======  ======  ======================
FASTQ_	Seq   
BED     Mapped  Skip mapping
BAM     Mapped  Skip mapping
======  ======  ======================

Cpipe may **not** support the following format in current version:

======  =======  ============================================
Format  Type     Solution
======  =======  ============================================
SRA   	Seq      Use `SRA Toolkit`_ to convert to FASTQ format
BED     Summit
BED     Peak
wig     Profile
Bigwig  Profile
======  =======  ============================================


Data preprocession
==================

Convert the raw sequencing data into intervals and profiles.

* Use Bowtie for tag alignment (mapping)
* Use MACS2 for peak calling


Correlation
===========

Focus on the visulization of similarity between replicates.

* Draw the venn diagram for peaks if there're less than 3 replicates (treatment or control)


Association Study
=================

Focus on association between intervals (result of peak calling) and traits like genome annotation.

* CEAS: Annotate the given intervals and scores with genome features 
* Conservation Plot: Calculates the PhastCons scores in several intervals sets

.. GO analysis
.. -----------

..   extract all the genes upstream or downstream the predicting peaks for functional clustering or annotation.


Motif
=====

  analysis the motif of the binding sites.

Quality control
===============
Based on Chip-seq pipeline and Cistrome DC database, QC program will generate a comprehensive quality control report about a particular dataset as well as the relative result compared to the whole DC database.

* Basic information: Species, Cell Type, Tissue Origin, Cell line, Factor, Experiment, Platform,  Treatment and Control. 
* Reads Genomic Mapping QC measurement: QC of raw sequence data with FastQC, FastQC score distribution, Basic mapping QC statistics, Mappable reads ratio, Mappable Redundant rate.
* Peak calling QC measurement: Peak calling summary, High confident Peak, Peaks overlapped with DHS(Dnase Hypersensitivity sites), Velcro ratio(human only), Profile correlation within union peak regions, Peaks overlap between Replicates.
* Functional Genomic QC measurement: Peak Height distribution, Meta Gene distribution, Peak conservation score, Motif QCmeasurement analysis.




.. _FASTQ: http://en.wikipedia.org/wiki/FASTQ_format


.. _SRA Toolkit: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software
