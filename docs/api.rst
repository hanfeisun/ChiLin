API Documentation
=================


Controllers of QC
------------------

.. automodule:: chilin.qc

.. autoclass:: chilin.qc.QC_Controller
    :members: runcmd, if_runcmd, _check, _render


.. autoclass:: chilin.qc.RawQC
    :members:  _infile_parse, _fastqc_info, run

.. autoclass:: chilin.qc.MappingQC
    :members: _basic_mapping_statistics_info, _mappable_ratio_info, _redundant_ratio_info, run

.. autoclass:: chilin.qc.PeakcallingQC
    :members: _peak_summary_info, _high_confidentPeaks_info, _velcro_ratio_info, _DHS_ratio_info, _replicate_info, run

.. autoclass:: chilin.qc.AnnotationQC
    :members: _ceas_info, DictToList, _distance, _conservation_info, get_seqpos, motif_info, run
        
.. autoclass:: chilin.qc.SummaryQC
    :members: run, packfile
    
DC class design instructions
--------------------------------

.. automodule:: chilin.dc

.. autoclass:: chilin.dc.LogWriter
     :members: record

.. autoclass::  chilin.dc.PipePreparation
     :members: get_config, get_rule, get_tex, get_summary, func_log, _read_rule, check_conf,

.. autoclass:: chilin.dc.PipeController
     :members: run_cmd, ifnot_runcmd, smart_run, cp, _render

.. autoclass:: chilin.dc.PipeGroom
     :members: run

.. autoclass:: chilin.dc.PipeBowtie
     :members: _sam2bam, _format, _extract, process

.. autoclass:: chilin.dc.PipeMACS2
     :members: _format, extract, run
 
.. autoclass:: chilin.dc.PipeVennCor
     :members: _format, extract, run

.. autoclass:: chilin.dc.PipeCEAS
     :members: _format, extract, run
 
.. autoclass:: chilin.dc.PipeConserv
     :members: _format, extract, _format, run

.. autoclass:: chilin.dc.PipeMotif
     :members: _format, extract, run

