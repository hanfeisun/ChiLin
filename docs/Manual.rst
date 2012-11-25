============
Manual
============
This manual is intended for two kind of users group.

**Biologists Users**: in order to know where they should go to download these
data and tools,what each output file represents

**Developers**: in order to make sure they're using the right format
of data and right version of tool to test, they have a consistent
naming convention and can find each file easily.

Through the pipeline, several temporary files will be generated, some of them are only used for settings
and transitions, others for continuing the next step, the rest for publishing and interpreting a biological
story. Below is three sections of tables for universal name rules.

.. note::
   %(DatasetID)s denotes *id* you input in the *basis* section in
   ChiLin.conf, %(treat_rep)s for treat number, %(control_rep)s for
   control rep. see :ref:`Get Started`
   use clear term to fill in the id
   example: use factor name plus your favorate number to replace the DatasetID below.
   if data is published, GSEID is recommended.

.. _Manual:

Raw Data
========

supported format
----------------------

single end *fastq* files and *absolid colored fastq* files are the raw reads files supported,
and users could use the *bam* files that has been already mapped back to the genomes.

ChiLin have been supporting the following format as input:

.. _supported formats<raw data>:
.. _Raw Data:

======  ======  ==========================================
Format  Type    instruction
======  ======  ==========================================
FASTQ_   Seq     single-end fastq or absolid colored fastq
BAM     Mapped  Skip mapping
======  ======  ==========================================

ChiLin do **not** support the following format in current version:

======  =======  ===============================================
Format  Type     Solution
======  =======  ===============================================
SRA     Seq      Use `SRA Toolkit`_ to convert to FASTQ format
BED     Summit
BED     Peak     Could be converted to bam files using `bedToBam`
wig     Profile  Use wiggleToBigwiggle to convert to Bigwiggle
Bigwig  Profile  Compressed Bigwiggle_
======  =======  ===============================================

Quality Control
---------------------
This part only use the tools for raw fastq quality control

1. tools involved
   We integrated Babraham's FastQC to assess the raw FastQ format
   files and extracted sequence quality scores from their summary,

2. Historic Data to compare
   new data will be stored for further comparing with all Cistrome_
   DC project historic data, which is collected by all DC team members
   and save into sqlite3 database format for indexing.

Output of Raw
------------------
1. temporary files

.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 30, 15
   :delim: ;

   FastQC summary ; %{DatasetID}s_fasctqc_summary.txt ; FastQC ; ref:`FastQC`
   FastQC treat data; %(DatasetID)s_rep%(treat_rep)s_treat_fastqc; FastQC_
   FastQC control data; %(DatasetID)s_rep%(control_rep)s_control_fastqc; FastQC_
   FastQC score R code; %(DatasetID)s_fastqc_score_distribution.r; R_ FastQC_
   FastQC score pdf ; %(DatasetID)s_fastqc_score_distribution.pdf; R_

2. No final result for package by default setting in this step

Reads Mapping
==================

Raw reads mapping is the first step for analyzing ChIP-seq data, which
is very important for following analysis.

Data analysis
---------------------
Modern high throughput sequencers can generate tens of mil-
lions of sequences in a single run. Before analysing this sequence
to draw biological conclusions you should always perform some
simple quality control checks to ensure that the raw data looks
good and there are no problems or biases in your data which
may a detect how you can usefully use it.

Here, we have chosen the bowtie for mapping raw reads data with
standard parameters.
below is the example command line we set in python script for `hg19`
::

  > /usr/local/bin/bowtie -p 1 -S -m 1 /pathtoIndex/hg19 /pathto/treat /outputdirectory/treat.sam

Quality Control
--------------------
Built-in tools would extract quality control preparation information
from standard output of Bowtie and sam files to do the following
description statistics.

Output of Mapping
-----------------
1. temporary files
.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 30, 15
   :delim: ;

   Bowie treat files ; %(DatasetID)s_treat_rep%{treat_rep}s.sam ; :ref:`Bowtie`
   Bowtie control files ; %(DatasetID)s_control_rep%{control_rep}s.sam ; :ref:`Bowtie`
   Bowtie temporary summary ; bowtie.tmp ; :ref: `Bowtie`

2. output files
Bowtie result sam files would be converted to bam binary format for
minimizing the file sizes through samtools::

    samtools view -bt chrom_len sam bam

.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 30, 15
   :delim; ;

   Bowtie treat bam ; %{DatasetID}s_%{control_rep}s.bam  ; :ref:`samtools`
   Bowtie control bam ; %{DatasetID}s_%{treat_rep}s.bam ; :ref:`samtools`

Groom
=====
This part is designed for users who don't have raw reads fastq files,
but have bam files instead. ChiLin helps to convert bam files to
fastq files for further processing all pipeline.
The only different for Usage is to input bam suffix files in the
*ChiLin.conf*

The convert tool used here is bedtools bamToFastq::

   bamToFastq -i x.bam -fq test.fq

Peak Calling
=============

Data Analysis
--------------------
We do the peak calling analysis by MACS2,
we set the parameter to meet the requirement for non redundant tags
for further analysis. Cutoff of false discovery rate(fdr) is set to
0.01, only keep one tag for duplicate tags to remove possible bias.
In default setting, macs2 would not build model.
The standard command line involves here::

     macs2 callpeak -B -q 0.01 --keep-dup 1 --shiftsize=73 --nomodel  -t /pathto/treat.bam  -c /pathto/control.bam -n macsname

Quality Control
-------------------
1. Before peaks calling, bams files would be sent for MACS2_ subparser
*filterdup* for statistical analysis on mapped reads non redundant
rate, the higher the measurement is, the better the data quality are.
2. After peaks calling, There are three measurement involves here,
total peaks count, confident peaks count , and shift size(optionally,
used when model=Yes, see in :ref:`advancedconf`)

Output files
------------
1. temporary files

.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 30, 15
   :delim: ;

   separate treat bedGraphfile ; %{DatasetID}s_treat_rep%{treat_rep}s.bdg ; :ref:`MACS2<MACS2>`
   separate control MACS bedGraph file ; %{DatasetID}s_control_rep%{control_rep}s.bdg ; :ref:`MACS2<MACS2>`
   Overall MACS bedGraph file ; %{DatasetID}s_treat.bdg ; :ref:`MACS2<MACS2>`
   bedGraph temporary file(remove exceptions) ; %{DatasetID}s_treat.bdg.tmp ;  :ref:`MACS2<MACS2>`
   sortedbed(For get top peaks) ; %(DatasetID)s_sorted.bed ; Linux sort
   top 1000 peaks(for latter MDSeqpos) ; %{DatasetID}s_top1000_summits.bed ; :ref:`MACS<MACS2>`
   MACS encode Peak(macs2 output) ; %{DatasetID}s_treat_rep%{treat_rep}s_peaks.encodePeak ; :ref:`MACS<MACS2>` 
   treatrep_pq_table ; %(DatasetID)s_rep%(treat_rep)s_pq_table.txt ; :ref:`MACS<MACS2>`
   pq_table ; %(DatasetID)s_pq_table.txt ; :ref:`MACS<MACS2>`
   treat_rep%(DatasetID)s_rep%(treat_rep)s_treat_pvalue.bdg; :ref:`MACS<MACS2>`
   treat_pvalue ; %(DatasetID)s_treat_pvalue.bdg; :ref:`MACS<MACS2>`
   treatrep_qvalue ; %(DatasetID)s_rep%(treat_rep)s_treat_qvalue.bdg ; :ref:`MACS<MACS2>`
   lambda_bdg ; %(DatasetID)s_rep%(control_rep)s_control_lambda.bdg ; :ref:`MACS<MACS2>`

2. final results

.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 25, 20
   :delim: ;

   treatreppeaks ; %(DatasetID)s_rep%(treat_rep)s_peaks.bed ; :ref:`MACS<MACS2>`
   treatpeaks ; %(DatasetID)s_peaks.bed ; :ref:`MACS<MACS2>`
   treatrepbw ; %(DatasetID)s_treat%(treat_rep)s.bw ; :ref:`MACS<MACS2>`
   treatbw ; %(DatasetID)s_treat.bw ; :ref:`MACS<MACS2>`
   controlrepbw ; %(DatasetID)s_rep%(treat_rep)s_control ; :ref:`MACS<MACS2>`
   controlbw ; %(DatasetID)s_control.bw ; :ref:`MACS<MACS2>`
   peaksrepxls ; %(DatasetID)s_rep%(treat_rep)s_peaks.xls ; :ref:`MACS<MACS2>`
   peaksxls ; %(DatasetID)s_peaks.xls ; :ref:`MACS<MACS2>`
   summitsrep ; %(DatasetID)s_rep%(treat_rep)s_summits.bed ; :ref:`MACS<MACS2>`
   summits ; %(DatasetID)s_summits.bed ; :ref:`MACS<MACS2>`

Replicates analysis
===================

Data analysis
-------------
Focus on the visulization of similarity between replicates.
* Draw the venn diagram for peaks if there're less than 3 replicates (treatment or control)
* plot the Correlation score for whole genome region average peaks score

Quality Control
---------------

The R code is searched by regular expression to get the needed part
for generating :ref:`QC report`

.. csv-table::
   :header: "Folder", "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 10
   :delim: ;
   
   correlation plot code ; %{DatasetID}s_cor.R ; :ref:`Buit-in tools`
   DHS peaks intersection ; %{DatasetID}s_bedtools_dhs.txt ; :ref:`BEDtools`
   overlap with velcro region ; %{DatasetID}s_bedtools_velcro.txt ; :ref:`BEDtools`
   peaks overlapped ; %{DatasetID}s_overlapped_bed ; :ref:`bedtools`
   AnnotationQC ; %{DatasetID}s_Metagene_distribution.pdf ; R
   AnnotationQC ; %{DatasetID}s_peak_height_distribution.pdf ; R

Meta genomics Study
=================

Focus on association between intervals (result of peak calling) and traits like genome annotation.

* CEAS: Annotate the given intervals and scores with genome features
* Conservation Plot: Calculates the PhastCons scores in several intervals sets

output files
------------
CEAS part 

.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 30, 15
   :delim: ;

    CEAS script ; %(DatasetID)s_ceas_CI.R ; :ref:`CEAS`
    CEAS script ; %(DatasetID)s_ceas_CI.pdf; :ref:`CEAS`
    CEAS xls ; %(DatasetID)s_ceas.xls; :ref:`CEAS`
    CEAS R script ; %(DatasetID)s_ceas.R; :ref:`CEAS`
    CEAS result pdf ; %(DatasetID)s_ceas.pdf

Conservation analysis  part

.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 30, 15
   :delim: ;

    conservtopsummits ; %(DatasetID)s_top3000summits.bed ; :ref:`built-in tools`
    conservR ; %(DatasetID)s_conserv.R ; :ref:`built-in tools`
    conservpng ; %(DatasetID)s_conserv.png ; :ref:`built-in tools`

Motif
=====
Here, we use the powerful combination of denovo motif finding
algorithm, MDscan, and database-based search algorithm, Seqpos for
motif analysis.

output files
------------
.. csv-table::
   :header: "Content", "File Name", "Tool used"
   :widths: 20, 30, 15
   :delim: ;

   summitspeaks1000 ; %(DatasetID)s_summits_p1000.bed; linux tools
   bgfreq ; %(DatasetID)s_bgfreq ; :ref:`MDSeqpos`
   seqpos ; %(DatasetID)s_seqpos.zip ; :ref: `MDSeqpos`

Other analysis type
===================

GO analysis
-----------

  extract all the genes upstream or downstream the predicting peaks for functional clustering or annotation.

Cistrome Radar/ Finder
----------------------
  You could check your top rated peaks in the Cistrome Radar and
  Finder to find interesting associated results, 

Summary Report
===================

Data Analysis Summary text
---------------------------

.. csv-table::
   :header: "Folder", "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 15
    root directory ; %{DatasetID}slog ; log; class Log

Quality report Instruction
--------------------------

An example QC report is here QC_.

.. Note:: 
   Output Format is optional(default PDF)
   Below is output in the root directory, that is the folder named after ${DatasetID}

Provide the overall report of the whole pipeline for viewing general result.


.. csv-table::
   :header: "Folder", "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 15
   :delim: ;
   root directory ; %{DatasetID}s_GSMID_QC.pdf ; All quality control measurements ; Main program

.. _PDF report:

.. _QC report:

Based on Chip-seq pipeline and Cistrome DC database, QC program will
generate a comprehensive quality control report about a particular
dataset as well as the relative result compared to the whole DC
database.
* QC report summary information ::
     Give an overview of all the measurement pass or fail information

* Basic information: Species, Cell Type, Tissue Origin, Cell line, Factor, Experiment, Platform,  Treatment and Control.
* Reads Genomic Mapping QC measurement: QC of raw sequence data with FastQC, FastQC score distribution, Basic mapping QC statistics, Mappable reads ratio, Mappable Redundant rate.
* Peak calling QC measurement: Peak calling summary, High confident Peak, Peaks overlapped with DHS(Dnase Hypersensitivity sites), Velcro ratio(human only), Profile correlation within union peak regions, Peaks overlap between Replicates.
* Functional Genomic QC measurement: Peak Height distribution, Meta Gene distribution, Peak conservation score, Motif QCmeasurement analysis.

.. _R: http://www.r-project.org/
.. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. _pdftex site: http://www.tug.org/applications/pdftex/
.. _samtools: samtools.sourceforge.net/SAM1.pdf
.. _Bigwiggle: http://genome.ucsc.edu/goldenPath/help/bigWig.html
.. _FASTQ: http://en.wikipedia.org/wiki/FASTQ_format
.. _SRA Toolkit: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software
.. _Processed Data:
.. _Cistrome: http://Cistrome.org
.. _MACS site: https://github.com/taoliu/MACS
.. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. _MDSeqPos site: https://bitbucket.org/cistrome/cistrome-applications-harvard/src/c477732c5c88/mdseqpos
.. _bedtools site: http://code.google.com/p/bedtools/
.. _SAMtools site: http://samtools.sourceforge.net/
.. _Bowtie site: http://bowtie-bio.sourceforge.net/index.shtml
.. _UCSC utilities: http://hgdownload.cse.ucsc.edu/admin/exe/
.. _UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables
.. _Cistrome site: http://cistrome.org/~hanfei
.. _FastQC site: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _pdfTex site: http://www.tug.org/applications/pdftex/ 
.. _QC: http://compbio.tongji.edu.cn/~meisl/document/7119_example.pdf
.. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. _pdftex site: http://www.tug.org/applications/pdftex/
