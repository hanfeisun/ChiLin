========================
Output & Temporary files
========================

This page is intended for:

**Developers**: in order to make sure they have a consistent naming convention and can find each file easily

**Users**: in order to know what each output file represents

Through the pipeline, several temporary files will be generated, some of them are only used for settings and transitions, others for continuing the next step, the rest for publishing and interpreting a biological story.Below is three sections of tables for universal name rules.

.. note::
     use clear term to replace the ${DatasetID}

     example: use factor name plus your favorate number to replace the DatasetID below.
     if data is published, GSEID is recommended.

Notation
========
All the program operation will be under the ${DatasetID} or 
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
   :header: "FolderName" "FileName", "Content", "Tool used"
   :widths: 20, 25, 20, 10
   :delim: ;
   
   ${DatasetID}_Bowtietmp ; ${DatasetID}_treat_rep${treat_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${DatasetID}_Bowtietmp ; ${DatasetID}_control_rep${control_rep}.sam ; mapping result ; :ref:`Bowtie`
   ${DatasetID}_Bowtietmp ; ${DatasetID}_bowtie_sh.txt ; bowtie shell summary : :ref: `Bowtie`
   ${DatasetID}_BEDtoolstmp ; ${DatasetID}_bedtools_dhs.txt ; DHS peaks intersection ; :ref:`BEDtools`
   ${DatasetID}_BEDtoolstmp ; ${DatasetID}_bedtools_velcro.txt ; overlap with velcro region; :ref:`BEDtools`
   ${DatasetID}_BEDtoolstmp ; ${DatasetID}_overlapped_bed ; peaks overlapped ; :ref:`bedtools`
   ${DatasetID}_MACStmp ; ${DatasetID}_control_rep${control_rep}.bdg ; separate control MACS bedGraph file; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat_rep${treat_rep}.bdg ; separate treat bedGraphfile ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat.bdg ; Overall MACS bedGraph file; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_treat.bdg.tmp ; bedGraph temporary file ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_rep${treat_rep}_treat.bdg ; separate treat bedGraph ; :ref:`MACS2<MACS2>`
   ${DatasetID}_MACStmp ; ${DatasetID}_${treat_rep}_peaks.encodePeak; MACS encode Peak ; :ref:`MACS<MACS2>`
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
   :header: "Folder" "File Name", "Content", "Tool used"
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
   ${DatasetID}_RawQCresult ; ${DatasetID}_fastqc_score_distribution.pdf ; Raw data QC ; R
   ${DatasetID}_PeakCallingQCresult ; ${DatasetID}_DHS_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_PeakCallingQCresult ; ${DatasetID}_velcro_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_PeakCallingQCresult ; ${DatasetID}_peak_ratio.pdf ; Peak calling QC ; R
   ${DatasetID}_QCresult ; ${DatasetID}_QC.tex ; QC report code ; pdftex_
   ${DatasetID}_QCresult ; ${DatasetID}_QC.tex ; QC report code ; :ref:`pdftex`
   root directory ; ${DatasetID}_summary.txt ; Data analysis summary ; : ref : `Built-in tools<Built-in tools>`

.. _PDF report:

Final PDF Report
================

Provide the overall report of the whole pipeline for viewing general result.
.. Note:: 
   Output Format is optional(default PDF)
   Below is output in the root directory, that is the folder named after ${DatasetID}

.. csv-table::
   :header: "Foler" "File Name", "Content", "Tool used"
   :widths: 20, 25, 20, 15
   :delim: ;

   root directory ; ${DatasetID}_ceas_combined.pdf  ; Cistron annotation ;  CEAS
   root directory ; ${DatasetID}_GSMID_QC.pdf ; All quality control measurements ; Main program
