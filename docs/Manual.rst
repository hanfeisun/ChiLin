============
Manual
============

This page is intended for:
**Users**: in order to know where they should go to download these
  data and tools,what each output file represents

**Developers**: in order to make sure they're using the right format
of data and right version of tool to test, they have a consistent
naming convention and can find each file easily

`Raw Data`
=============

supported format
----------------------

single end *fastq* files and *absolid colored fastq* files are the raw reads files supported,
and users could use the *bam* files that has been already mapped back to the genomes.

ChiLin have been supporting the following format as input:

======  ======  ==========================================
Format  Type    instruction
======  ======  ==========================================
FASTQ   Seq     single-end fastq or absolid colored fastq
BAM     Mapped  Skip mapping
======  ======  ==========================================

ChiLin may **not** support the following format in current version:

======  =======  ===============================================
Format  Type     Solution
======  =======  ===============================================
SRA     Seq      Use `SRA Toolkit`_ to convert to FASTQ format
BED     Summit   
BED     Peak     Could be converted to bam files using `bedToBam`
wig     Profile
Bigwig  Profile
======  =======  ===============================================

Quality Control
---------------------
1. tools involved

2. Historic Data to compare


Reads Mapping 
==================

Data analysis
---------------------

Quality Control
--------------------


Peak Calling
=============

Data Analysis
--------------------

Quality Control
-------------------


Summary Report
===================

Data Analysis Summary text
---------------------------


Quality Control Report
--------------------------

Provide the overall report of the whole pipeline for viewing general result.

.. Note:: 
   Output Format is optional(default PDF)
   Below is output in the root directory, that is the folder named after ${DatasetID}

.. csv-table::
   :header: "Folder", "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 15
   :delim: ;

   root directory ; ${DatasetID}_ceas_combined.pdf  ; Cistron annotation ;  CEAS
   root directory ; ${DatasetID}_GSMID_QC.pdf ; All quality control measurements ; Main program
.. _PDF report:


.. _Raw Data:



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

.. ====
.. Data
.. ====

.. Built-in Data
.. -------------

.. The Cpipe package includes all the build-in data for hg19 and mm9. For other species, you may need to download these data from data source or custom it yourself.
 
.. ============================   ============  =====================  =========  
.. Data Name                       Used by       Data Source           Format     
.. ============================   ============  =====================  =========  
.. Chromesome length              samtools      `UCSC table browser`_  2-column   
.. Chromesome length              CEAS          --                     --
.. Genome backgroud annotation    CEAS          `CEAS site`_           sqlite3
.. DHS region                     bedtools      Custom                 BED
.. Velcro region                  bedtools	     Custom                 BED
.. Motif database                 MDSeqPos      `MDSeqPos site`_       xml
.. FastQC result database         QCreport      Custom                 bed
.. Data summary database          QCreport      Custom                 bed
.. ============================   ============  =====================  =========


.. .. _External Data:

.. External Data
.. -------------

.. Some data are too large to be included by the pipeline package, so you need to download these data from data source.

.. ============================   =================  =====================  =========  
.. Data Name                       Used by           Data Source            Format     
.. ============================   =================  =====================  =========  
.. Bowtie pre-built index         Bowtie             `Bowtie site`_         ebwt
.. Conservation profile           Conservation Plot  `Cistrome site`_       Bigwig
.. ============================   =================  =====================  =========  

.. =====
.. Tools
.. =====

.. Built-in Tools
.. --------------

.. Built-in tools are the scripts that can be run from command-line independently when you have installed the Cpipe package.


.. .. _Built-in tools:

.. ============================   =====================  
.. Tool Name                      Modified from        
.. ============================   =====================  
.. liftover
.. Venn Diagram
.. Conservation Plot
.. Correlation plot               bigwig_correlation
.. bamtofastq
.. BedClip
.. wigTobigwiggle
.. RegPotential
.. sample_contamination
.. ============================   =====================  


.. .. _Bowtie:
.. .. _samtools:
.. .. _MACS2:
.. .. _MDSeqpos:
.. .. _BEDtools:
.. .. _External Tools:

.. External Tools
.. --------------


.. External Tools are the tools invoked by Cpipe by their path.

.. ============================   =====================  ==================    
.. Tool Name                      Download source         Version
.. ============================   =====================  ==================    
.. FastQC
.. R
.. Cython
.. MACS2                          `MACS site`_           2.0.10 20120605
.. CEAS                           `CEAS site`_           0.9.9.7
.. bedtools		       `bedtools site`_	      v2.16.2
.. pybedtools
.. samtools		       `SAMtools site`_	      0.1.17
.. Bowtie                         `Bowtie site`_         0.12.8
.. bedGraphToBigWig	       `UCSC utilities`_      v4
.. FastQC                         `FastQC site`_         v0.10.1
.. pdfTeX                         `pdfTex site`_         v1.40.10
.. IGV
.. ============================   =====================  ==================    


.. .. _MACS site: https://github.com/taoliu/MACS
.. .. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. .. _MDSeqPos site: https://bitbucket.org/cistrome/cistrome-applications-harvard/src/c477732c5c88/mdseqpos
.. .. _bedtools site: http://code.google.com/p/bedtools/
.. .. _SAMtools site: http://samtools.sourceforge.net/
.. .. _Bowtie site: http://bowtie-bio.sourceforge.net/index.shtml
.. .. _UCSC utilities: http://hgdownload.cse.ucsc.edu/admin/exe/
.. .. _UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables
.. .. _Cistrome site: http://cistrome.org/~hanfei
.. .. _FastQC site: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. .. _pdfTex site: http://www.tug.org/applications/pdftex/ 

.. ========
.. Workflow
.. ========

.. .. digraph:: foo

..     rankdir=TB
..     size="15,15"
..     edge[arrowhead=open]

..     start[shape=circle, label="", style=filled]
..     end[shape=doublecircle, label="", style=filled]

..     readconf[shape=box,style=rounded, label="class Check"]
..     bowtie[shape=box,style=rounded, label="Run Bowtie"]
..     rawQC[shape=box,style=rounded, label="Run RawQC"]
..     mappingQC[shape=box,style=rounded, label="Run MappingQC"]
..     macs2[shape=box,style=rounded, label="Run MACS2"]
..     peakcallingQC[shape=box,style=rounded, label="Run PeakcallingQC"]
..     ceas_seqpos[shape=box,style=rounded, label="Run CEAS/Seqpos"]
..     venn[shape=box,style=rounded, label="class Replicates, Draw VennDiagram and Correlation plot"]
..     conservation[shape=box,style=rounded, label="Draw ConservationPlot"]
..     annotationQC[shape=box,style=rounded, label="Run AnnotationQC"]

    
..     ifmapped[shape=diamond, label="Mapped?"]
..     ifrep[shape=diamond, label="Replicate?"]
    
..     start -> readconf
..     readconf -> rawQC
..     rawQC -> ifmapped[headport=n, color="grey"]
..     ifmapped -> mappingQC[label="[Yes]" tailport=s]
..     ifmapped -> bowtie[taillabel="[No]" tailport=e]
..     bowtie -> mappingQC
..     mappingQC -> macs2[color="grey"]
..     macs2 -> ifrep
..     peakcallingQC -> ceas_seqpos[color="grey"]
..     ifrep -> venn[label="[Yes]" tailport=s]
..     ifrep -> conservation[label="[No]" tailport=e]
..     venn -> conservation
..     conservation -> peakcallingQC
..     ceas_seqpos -> annotationQC
..     annotationQC -> end[taillabel="Output Report"]


========================
Output & Temporary files
========================

This page is intended for:


Through the pipeline, several temporary files will be generated, some of them are only used for settings and transitions, others for continuing the next step, the rest for publishing and interpreting a biological story.Below is three sections of tables for universal name rules.

.. note::
     use clear term to replace the ${DatasetID}

     example: use factor name plus your favorate number to replace the DatasetID below.
     if data is published, GSEID is recommended.

Notation
========
All the program operation will be under the ${DatasetID}_folder.  

.. envvar:: ${DatasetID}

    The value of :ref:`dataset.id<dataset.id>` option in :envvar:`[meta]` section

.. envvar:: ${treat_rep}

    The suffix of :envvar:`treatment` option in :envvar:`[meta]` section


.. envvar:: ${control_rep}

    The suffix of :envvar:`control` option in :envvar:`[meta]` section

.. envvar:: ${config}

    The general configuration file for pipeline :envvar:`[meta]` section

.. envvar:: ${log}

    For write in all shell output and assessment during procedure, including time consumed :envvar: `[meta]`

Temporary files
===============

.. csv-table::
   :header: "FolderName", "FileName", "Content", "Tool used"
   :widths: 25, 25, 20, 10
   :delim: ;
   
   root directory ; ${DatasetID}log ; log; class Log
   ${DatasetID}_Bowtietmp ; ${DatasetID}_treat_rep${treat_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${DatasetID}_Bowtietmp ; ${DatasetID}_treat_rep${treat_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${DatasetID}_Bowtietmp ; ${DatasetID}_control_rep${control_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${DatasetID}_Bowtietmp ; ${DatasetID}_bowtie_sh.txt ; bowtie shell summary ; :ref: `Bowtie`
   ${DatasetID}_BEDtoolstmp ; ${DatasetID}_bedtools_dhs.txt ; DHS peaks intersection ; :ref:`BEDtools`
   ${DatasetID}_BEDtoolstmp ; ${DatasetID}_bedtools_velcro.txt ; overlap with velcro region; :ref:`BEDtools`
   ${DatasetID}_BEDtoolstmp ; ${DatasetID}_overlapped_bed ; peaks overlapped ; :ref:`bedtools`
   ${DatasetID}_MACStmp ; ${DatasetID}_control_rep${control_rep}.bdg ; separate control MACS bedGraph file; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat_rep${treat_rep}.bdg ; separate treat bedGraphfile ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat.bdg ; Overall MACS bedGraph file ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat.bdg.tmp ; bedGraph temporary file ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_treat.bdg ; separate treat bedGraph ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_${treat_rep}_peaks.encodePeak ; MACS encode Peak ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_pq_table.txt ; separate p q value  ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_pq_table.txt ; collective MACS2 p q value ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_control_lambda.bdg ; treat over control lambda; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_control.bdg ; treat over control ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_peaks.xls ; peaks calling list ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat_peaks.xls ; overall peak file ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_pq_table.txt ; peaks calling p q value ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_summits.bed ; peaks summits ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_treat_logLR.bdg ; log bedGraph ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat_logLR.bdg ; log bedGraph ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_treat_pvalue.bdg ; treat bedGraph pvalue ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat_pvalue.bdg ; treat overall p value ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_treat_qvalue.bdg ; treat bedGraph q value ;  :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_top1000_summits.bed ; top 1000 peaks ; :ref:`MACS<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_bgfreq ; MACS background frequence ; :ref:`MACS<MACS2>`
   ${DatasetID}_Cortmp ; ${DatasetID}_cor.R ; correlation plot code ; :ref:`Buit-in tools`
   ${DatasetID}_CEAStmp ; ${DatasetID}_ceaswithoutpeak.R ; CEAS ; R
   ${DatasetID}_CEAStmp ; ${DatasetID}_ceaswithpeak.R ; CEAS ; R
   ${DatasetID}_CEAStmp ; ${DatasetID}_ceaswithoutpeak.pdf ; CEAS ; R
   ${DatasetID}_CEAStmp ; ${DatasetID}_ceaswithpeak.pdf ; CEAS ; R
   ${DatasetID}_qctmp ; ${DatasetID}_fasctqc_summary.txt ; FastQC ; ref:`FastQC`
   ${DatasetID}_qctmp ; ${DatasetID}_Metagene_distribution.pdf ; AnnotationQC ; R
   ${DatasetID}_qctmp ; ${DatasetID}_peak_height_distribution.pdf ; AnnotationQC ; R


.. _Processed Data:

Output result
=============

.. csv-table:: 
   :header: "Folder", "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 10
   :delim: ;
   
   root directory ; ${DatasetID}log ; log; class Log
   ${DatasetID}_bowtieresult ; ${DatasetID}_${control_rep}.bam ; mapping result ; :ref:`samtools`
   ${DatasetID}_bowtieresult ; ${DatasetID}_${treat_rep}.bam ; mapping result; 
   ${DatasetID}_MACSresult ; ${DatasetID}_${treat_rep}_peaks.bed ;Peak calling ; :ref:`MACS2<MACS2>`      
   ${DatasetID}_corresult ; ${DatasetID}_cor.R ; correlation plot code ; :ref:`Built-in tools<Built-in tools>`
   ${DatasetID}_corresult ; ${DatasetID}_cor.pdf ; correlation plot pdf ; :ref:`Built-in tools<Built-in tools>`
   ${DatasetID}_Motifresult ; ${DatasetID}_seqpos.zip ; Motif analysis ; :ref:`MDSeqpos<MDSeqpos>`
   ${DatasetID}_CEASresult ;${DatasetID}_ceas.xls ; CEAS ; CEAS_
   ${DatasetID}_conservresult ; ${DatasetID}_conserv.png ; Phascon score plot ; :ref:`Built-in tools<Built-in tools>`
   ${DatasetID}_conservresult ; ${DatasetID}_conserv.R ; Phascon score ; :ref:`Built-in tools<Built-in tools>`
   ${DatasetID}_MappingQCresult ; ${DatasetID}_redundant_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_MappingQCresult ; ${DatasetID}_mappable_ratio.pdf ; Mapping QC result ; R
   ${DatasetID}_QCresult ; ${DatasetID}_fastqc_score_distribution.pdf ; Raw data QC ; R
   ${DatasetID}_QCresult ; ${DatasetID}_fastqc_summary.txt ; Raw data QC ; R
   ${DatasetID}_QCresult ; ${DatasetID}_DHS_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_QCresult ; ${DatasetID}_velcro_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_QCresult ; ${DatasetID}_peak_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_QCresult ; ${DatasetID}_QC.tex ; QC report code ; pdftex_
   ${DatasetID}_QCresult ; ${DatasetID}_QC.pdf ; QC report ; :ref:`pdftex`
   root directory ; ${DatasetID}_summary.txt ; Data analysis summary ; : ref : `Built-in tools<Built-in tools>`




.. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. _pdftex site: http://www.tug.org/applications/pdftex/
