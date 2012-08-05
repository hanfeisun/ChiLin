
API Documentation
=================


Controllers of QC
------------------

.. automodule:: chilin.qc

.. autoclass:: chilin.qc.QC_Controller
    :members: run, _check, _render


.. autoclass:: chilin.qc.RawQC
	:members:  _infile_parse, _fastqc_info, run

.. autoclass:: chilin.qc.MappingQC
	:members: _basic_mapping_statistics_info, _mappable_ratio_info, _redundant_ratio_info, run

.. autoclass:: chilin.qc.PeakcallingQC
	:members: _peak_summary_info, _high_confidentPeaks_info, _velcro_ratio_info, _DHS_ratio_info, _replicate_info, run

.. autoclass:: chilin.qc.AnnotationQC
	:members: _ceas_info, _distance, _conservation_info, _DictToList, _motif_info, run
.. autoclass:: chilin.qc.SummaryQC
	:members: run, packfile
    
DC class design instructions
--------------------------------

.. automodule:: chilin.dc

.. autoclass::  chilin.dc.PipePreparation
     :members: readconf, checkconf

.. autoclass:: chilin.dc.PathFinder
    :members: readconf, parseconfrep

.. autoclass:: chilin.chilin.LogWriter
     :members: record

.. autoclass:: chilin.dc.PipeController
     :members: run, partition, render

.. autoclass:: chilin.dc.PipeBowtie
     :members: _format, extract, process
     
.. autoclass:: chilin.dc.PipeMACS2
     :members: _format, _median, extract, process

 
.. autoclass:: chilin.dc.PipeVennCor
     :members: _format, extract, process

.. autoclass:: chilin.dc.PipeCEAS
     :members: _format, extract, process
 
.. autoclass:: chilin.dc.PipeConserv
     :members: _format, extract, process

.. autoclass:: chilin.dc.PipeMotif  
     :members: _format, extract, process

