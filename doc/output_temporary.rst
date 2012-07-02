========================
Output & Temporary files
========================

This page is intended for:

**Developers**: in order to make sure they have a consistent naming convention and can find each file easily

**Users**: in order to know what each output file represents

Through the pipeline, several temporary files will be generated, some of them are only used for settings and transitions, others for continuing the next step, the rest for publishing and interpreting a biological story.Below is three sections of tables for universal name rules.

.. note::
     For data which has not been published on GEO
     use factor name plus your favorate number to replace the GSMID below.

Notation
========

.. envvar:: ${DatasetID}

    The value of :ref:`dataset.id<dataset.id>` option in :envvar:`[meta]` section

.. envvar:: ${treat_rep}

    The suffix of :envvar:`treatment` option in :envvar:`[meta]` section


.. envvar:: ${control_rep}

    The suffix of :envvar:`control` option in :envvar:`[meta]` section    

Temporary files
===============

.. csv-table:: 
   :header: "Name", "Content", "Tool used"
   :widths: 20, 20, 10
   :delim: ;
   
   ${DatasetID}_treat_rep${treat_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${DatasetID}_control_rep${control_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${DatasetID}.conf ; configuration ; Main program
   ${DatasetID}_bedtools_dhs.txt ; DHS peaks intersection ; :ref:`BEDtools`
   ${DatasetID}_cor.R ; correlation plot code ; :ref:`Buit-in tools`
   ${DatasetID}_seqpos.zip ; Motif analysis ; :ref:`MDSeqpos`
   ${DatasetID}_QC.tex ; QC report code ; :ref:`pdftex`
   ${DatasetID}_mappable_ratio.pdf ; Mapping QC result ; R
   ${DatasetID}_fastqc_score_distribution.pdf ; Raw data QC ; R
   ${DatasetID}_peak_distribution.pdf ; Peak calling QC ; R
   ${DatasetID}_velcro_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_peak_overlap_DHS.pdf ; Peak calling QC ; R

.. _Processed Data:

Output result
=============

.. csv-table:: 
   :header: "Name", "Content", "Tool used"
   :widths: 20, 20, 10
   :delim: ;
   
   Folder ; containing all results ; Main  Program     
   GSMIDlog ; log;                                               
   GSMID_control_repnumber.bam ; mapping result ; :ref:`samtools`
   GSMID_treat_repnumber.bam ; mapping result; 
   GSMID_repnumber_peaks.bed ;Peak calling ; :ref:`MACS2<MACS2>`      
   GSMID_bedtools_dhs.txt ; DHS peaks intersection ; :ref:`BEDtools<BEDtools>`
   GSMID_cor.R ; correlation plot code ; :ref:`Built-in tools<Built-in tools>`
   GSMID_seqpos.zip ; Motif analysis ; :ref:`MDSeqpos<MDSeqpos>`
   GSMID_QC.tex ; QC report code ; pdftex_
   GSMID_ceas.xls ; CEAS ; CEAS_
   GSMID_conserv.png ; Phascon score plot ; :ref:`Built-in tools<Built-in tools>`
   GSMID_conserv.R ; Phascon score ; :ref:`Built-in tools<Built-in tools>`



.. _PDF report:

Final PDF Report
================

Provide the overall report of the whole pipeline for viewing general result.

+--------------------------+---------------------------------+--------------------+
| Name                     | Content                         |   Tool used        |
+==========================+=================================+====================+
| GSMID_ceas_combined.pdf  | Cistron annotation              |  CEAS              |
+--------------------------+---------------------------------+--------------------+
| GSMID_QC.pdf             | All quality control measurements| Main program       |
+--------------------------+---------------------------------+--------------------+
