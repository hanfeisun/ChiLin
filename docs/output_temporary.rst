========================
Output & Temporary files
========================

This page is intended for:

**Developers**: in order to make sure they have a consistent naming convention and can find each file easily

**users**: in order to know what each output file represents

Through the pipeline, several temporary files will be generated, some of them are only used for settings and transitions, others for continuing the next step, the rest for publishing and interpreting a biological story.Below is three sections of tables for universal name rules.

.. note::
     use clear term to replace the ${dataset_id}

     example: use factor name plus your favorate number to replace the dataset_id below.
     if data is published, GSEID is recommended.

Notation
========
All the program operation will be under the ${dataset_id}_folder.

.. envvar:: ${dataset_id}

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
   
   root directory ; ${dataset_id}log ; log; class Log
   ${dataset_id}_Bowtietmp ; ${dataset_id}_treat_rep${treat_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${dataset_id}_Bowtietmp ; ${dataset_id}_treat_rep${treat_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${dataset_id}_Bowtietmp ; ${dataset_id}_control_rep${control_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${dataset_id}_Bowtietmp ; ${dataset_id}_bowtie_sh.txt ; bowtie shell summary ; :ref: `Bowtie`
   ${dataset_id}_BEDtoolstmp ; ${dataset_id}_bedtools_dhs.txt ; DHS peaks intersection ; :ref:`BEDtools`
   ${dataset_id}_BEDtoolstmp ; ${dataset_id}_bedtools_velcro.txt ; overlap with velcro region; :ref:`BEDtools`
   ${dataset_id}_BEDtoolstmp ; ${dataset_id}_overlapped_bed ; peaks overlapped ; :ref:`bedtools`
   ${dataset_id}_MACStmp ; ${dataset_id}_control_rep${control_rep}.bdg ; separate control MACS bedGraph file; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_treat_rep${treat_rep}.bdg ; separate treat bedGraphfile ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_treat.bdg ; Overall MACS bedGraph file ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_treat.bdg.tmp ; bedGraph temporary file ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_treat.bdg ; separate treat bedGraph ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_${treat_rep}_peaks.encodePeak ; MACS encode Peak ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_pq_table.txt ; separate p q value  ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_pq_table.txt ; collective MACS2 p q value ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_control_lambda.bdg ; treat over control lambda; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_control.bdg ; treat over control ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_peaks.xls ; peaks calling list ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_treat_peaks.xls ; overall peak file ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_pq_table.txt ; peaks calling p q value ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_summits.bed ; peaks summits ; :ref:`MACS2<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_treat_logLR.bdg ; log bedGraph ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_treat_logLR.bdg ; log bedGraph ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_treat_pvalue.bdg ; treat bedGraph pvalue ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_treat_pvalue.bdg ; treat overall p value ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_rep${treat_rep}_treat_qvalue.bdg ; treat bedGraph q value ;  :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_top1000_summits.bed ; top 1000 peaks ; :ref:`MACS<MACS2>`
   ${dataset_id}_MACStmp ; ${dataset_id}_bgfreq ; MACS background frequence ; :ref:`MACS<MACS2>`
   ${dataset_id}_Cortmp ; ${dataset_id}_cor.R ; correlation plot code ; :ref:`Buit-in tools`
   ${dataset_id}_CEAStmp ; ${dataset_id}_ceaswithoutpeak.R ; CEAS ; R
   ${dataset_id}_CEAStmp ; ${dataset_id}_ceaswithpeak.R ; CEAS ; R
   ${dataset_id}_CEAStmp ; ${dataset_id}_ceaswithoutpeak.pdf ; CEAS ; R
   ${dataset_id}_CEAStmp ; ${dataset_id}_ceaswithpeak.pdf ; CEAS ; R
   ${dataset_id}_qctmp ; ${dataset_id}_fasctqc_summary.txt ; FastQC ; ref:`FastQC`
   ${dataset_id}_qctmp ; ${dataset_id}_Metagene_distribution.pdf ; AnnotationQC ; R
   ${dataset_id}_qctmp ; ${dataset_id}_peak_height_distribution.pdf ; AnnotationQC ; R


.. _Processed Data:

Output result
=============

.. csv-table:: 
   :header: "Folder", "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 10
   :delim: ;
   
   root directory ; ${dataset_id}log ; log; class Log
   ${dataset_id}_bowtieresult ; ${dataset_id}_${control_rep}.bam ; mapping result ; :ref:`samtools`
   ${dataset_id}_bowtieresult ; ${dataset_id}_${treat_rep}.bam ; mapping result;
   ${dataset_id}_MACSresult ; ${dataset_id}_${treat_rep}_peaks.bed ;Peak calling ; :ref:`MACS2<MACS2>`
   ${dataset_id}_corresult ; ${dataset_id}_cor.R ; correlation plot code ; :ref:`Built-in tools<Built-in tools>`
   ${dataset_id}_corresult ; ${dataset_id}_cor.pdf ; correlation plot pdf ; :ref:`Built-in tools<Built-in tools>`
   ${dataset_id}_Motifresult ; ${dataset_id}_seqpos.zip ; Motif analysis ; :ref:`MDSeqpos<MDSeqpos>`
   ${dataset_id}_CEASresult ;${dataset_id}_ceas.xls ; CEAS ; CEAS_
   ${dataset_id}_conservresult ; ${dataset_id}_conserv.png ; Phascon score plot ; :ref:`Built-in tools<Built-in tools>`
   ${dataset_id}_conservresult ; ${dataset_id}_conserv.R ; Phascon score ; :ref:`Built-in tools<Built-in tools>`
   ${dataset_id}_MappingQCresult ; ${dataset_id}_redundant_ratio.pdf ; Peak calling QC ; R
   ${dataset_id}_MappingQCresult ; ${dataset_id}_mappable_ratio.pdf ; Mapping QC result ; R
   ${dataset_id}_QCresult ; ${dataset_id}_fastqc_score_distribution.pdf ; Raw data QC ; R
   ${dataset_id}_QCresult ; ${dataset_id}_fastqc_summary.txt ; Raw data QC ; R
   ${dataset_id}_QCresult ; ${dataset_id}_DHS_ratio.pdf ; Peak calling QC ; R
   ${dataset_id}_QCresult ; ${dataset_id}_velcro_ratio.pdf ; Peak calling QC ; R
   ${dataset_id}_QCresult ; ${dataset_id}_peak_ratio.pdf ; Peak calling QC ; R
   ${dataset_id}_QCresult ; ${dataset_id}_QC.tex ; QC report code ; pdftex_
   ${dataset_id}_QCresult ; ${dataset_id}_QC.pdf ; QC report ; :ref:`pdftex`
   root directory ; ${dataset_id}_summary.txt ; Data analysis summary ; : ref : `Built-in tools<Built-in tools>`

.. _PDF report:

Final PDF Report
================
Provide the overall report of the whole pipeline for viewing general result.

.. Note:: 
   Output Format is optional(default PDF)
   Below is output in the root directory, that is the folder named after ${dataset_id}

.. csv-table::
   :header: "Folder", "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 15
   :delim: ;

   root directory ; ${dataset_id}_ceas_combined.pdf  ; Cistron annotation ;  CEAS
   root directory ; ${dataset_id}_GSMID_QC.pdf ; All quality control measurements ; Main program

.. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. _pdftex site: http://www.tug.org/applications/pdftex/
