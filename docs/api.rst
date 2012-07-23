
API Documentation
=================

DC class design instructions
--------------------------------

.. automodule:: chilin.dc

.. autoclass::  chilin.dc.PipePreparation
     :members: readconf, checkconf

.. autoclass:: chilin.dc.PathFinder
    :members: bowtiefilepath, macs2filepath, venn_corfilepath, ceasfilepath, conservfilepath, qcfilepath

.. autoclass:: chilin.chilin.LogWriter
     :members: record

.. autoclass:: chilin.dc.PipeController
     :members: run, partition, render

.. autoclass:: chilin.dc.PipeBowtie
     :members: _format, _run, summary
     
.. autoclass:: chilin.dc.PipeMACS2
     :members: _format, _run, summary

 
.. autoclass:: chilin.dc.PipeVennCor
     :members: _format, _run, summary

.. autoclass:: chilin.dc.PipeCEAS
     :members: _format, _run, summary
 
.. autoclass:: chilin.dc.PipeConserv
     :members: _format, _run, summary

.. autoclass:: chilin.dc.PipeMotif  
     :members: _format, _run, summary


Controllers of QC
------------------

.. automodule:: chilin.qc

.. autoclass:: chilin.qc.QC_Controller
    :members: run, check, render


.. autoclass:: chilin.qc.RawQC
	:members: _basic_info, _fastqc_info, run, check

.. autoclass:: chilin.qc.MappingQC
	:members: _basic_mapping_statistics_info, _mappable_ratio_info, _redundant_ratio_info, run, check

.. autoclass:: chilin.qc.PeakcallingQC
	:members: _peak_summary_info, _velcro_ratio_info, _DHS_ratio_info, _replicate_info, run, check

.. autoclass:: chilin.qc.AnnotationQC
	:members: _ceas_info, _conservation_info, _motif_info, run, check
    
