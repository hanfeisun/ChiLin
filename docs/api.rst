
API Documentation
=================

DC class design instructions
--------------------------------

.. automodule:: chilin.dc

.. autoclass::  chilin.dc.Check
     :member: ReadConf, CheckConf, DependencyCheck

.. autoclass:: chilin.dc.TemplateParser
    :member: Loader

.. autoclass:: chilin.chilin.Log
     :member: _Timer, info

.. autoclass:: chilin.dc.DcController
     :member:_StepControl, run, partition, render

.. autoclass:: chilin.dc.Bowtie
     :member: _format, _run, summary
     
.. autoclass:: chilin.dc.MACS
     :member: _format, _run, summary

 
.. autoclass:: chilin.dc.Replicates
     :member: _format, _run, summary

.. autoclass:: chilin.dc.CEAS
     :member: _format, _run, summary
 
.. autoclass:: chilin.dc.Conserv
     :member: _format, _run, summary



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
    
