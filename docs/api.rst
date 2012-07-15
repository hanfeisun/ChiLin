
API Documentation
=================




Controllers of QC
------------------

.. automodule:: chilin.qc

.. .. automodule:: chilin.chilin

.. .. autoclass::  chilin.chilin.Check
..     :member: ReadConf, CheckConf, DependencyCheck

.. .. autoclass:: chilin.chilin.TemplateParser
..     :member: Loader, DcRender, QcRender

.. .. autoclass:: chilin.chilin.Log
..     :member: _Timer, info

.. .. autoclass:: chilin.chilin.DcController
..     :member:_StepControl, main
.. 
.. .. autoclass:: chilin.chilin.Bowtie
..     :member: _format, _run, summary
..     
.. .. autoclass:: chilin.chilin.MACS
..     :member: _format, _run, summary
.. 
.. .. autoclass:: chilin.chilin.Bowtie
..     :member: _format, _run, summary
.. 
.. .. autoclass:: chilin.chilin.Bowtie
..     :member: _format, _run, summary
.. 
.. .. autoclass:: chilin.chilin.Conserv
    :member: _format, _run, summary

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
    
