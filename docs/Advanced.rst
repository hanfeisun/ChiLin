==============================
Advanced sections for ChiLin
==============================
This sections for developers who are interested in our development
and want to participate
**Developers**: in order to be sure they are developing the right project that fulfills requirements provided in this document.
**Users**: in order to get familiar with the idea of the project and suggest other features that would make it even more functional. 

Dive into options
===================

Configuration template generator
---------------------------------
The template use jinja method of assigning methods, the template
system is extensible and flexible, you could customize your desired
species analysis dependent conf.
Just modify the ChiLin.conf under chilin/conf *{% set %}* part.

step control
---------------

If you don't want to go from beginning to the end of the pipeline, say, you just want the peaks calling results,
the choice is the *-d* option, you could end at wherever you want.

Debug mode
-------------
Here, we design options for debug model by using *--debug*
Through checking whether the output files are in the outputdirectory, the program will judge by itself
to continue or skip the steps.

only quality control
--------------------
-k options for only quality control
if you have run dc pipeline, or you just want to see the preprocessed data
quality, you can run qc only with *-k* option

Multithread
--------------
For Bowtie and Fastqc, we provide options *-p* to open inner multi pool
supporting much faster processing.

Hidden features
=================

Macs model or not in the conf files

if further analysis requirements is needed, you may want to ajust the
default options we set here.

.. _advancedconf:


Dive into conf
======================



for ceas
-------------------
For Dnase seq, the default value for peaks number is all because of
the unique data pattern. For Histone and TF, the default value is set
to 3000.

* bipromotor sizes
* promotor sizes
* rel_dist ::

for conservation
------------------
* w width of conservation plot region

for macs2
-----------------
* --shift-size --nomodel: optional

motif
-----
1. MDscan width is the speed limited step in the denovo motif
discovery.

conf as following ::
   [seqpos]
   ...
   mdscan_width = 200

2.

Configuration instructions
----------------------------

.. envvar:: [basis]

    Lists all the meta-data of current workflow.
    Consist of the following options:

    .. envvar:: id

        The name for the dataset, which will be the value of :envvar:`%{DatasetID}s`
        Limit: a string (1) consist of ``numbers``, ``alphabets`` or ``'_'`` (2) shorter than 20 characters

    .. envvar:: species
        The name of species, written to the QCreport and log
        Limit: a string (1) consist of ``numbers``, ``alphabets`` or ``'_'`` (2) shorter than 20 characters

    .. envvar:: factor

        The name of species, writen to DC summary and QCreport, log
        Limit: a string (1) come from GO standard term


    .. envvar:: treat

       The paths of treatment files
       Limit: absolute ``path`` of files in :ref:`supported formats<raw data>`

    .. envvar:: control

       The paths of treatment files
       Limit: absolute or relative ``path`` of files in :ref:`supported formats<raw data>`

Dive into rule
===============
We separate Name Rules for output and temporary files from analysis codes part for easier to maintain,
If you don't feel comfortable of our name ways, it's simple for you to adjust it.


Dive into background
====================

Built-in Data
-------------

The ChiLin package includes all the build-in data for hg19 and mm9. For other species, you may need to download these data from data source or custom it yourself.

============================   ============  =====================  =========  
Data Name                       Used by       Data Source           Format     
============================   ============  =====================  =========  
Chromesome length              samtools      `UCSC table browser`_  2-column   
Chromesome length              CEAS          --                     --
Genome backgroud annotation    CEAS          `CEAS site`_           sqlite3
DHS region                     bedtools      Custom                 BED
Velcro region                  bedtools	     Custom                 BED
Motif database                 MDSeqPos      `MDSeqPos site`_       xml
FastQC result database         QCreport      Custom                 bed
Data summary database          QCreport      Custom                 bed
============================   ============  =====================  =========


.. _External Data:

External Data
-------------

Some data are too large to be included by the pipeline package, so you need to download these data from data source.

============================   =================  =====================  =========  
Data Name                       Used by           Data Source            Format     
============================   =================  =====================  =========  
Bowtie pre-built index         Bowtie             `Bowtie site`_         ebwt
Conservation profile           Conservation Plot  `Cistrome site`_       Bigwig
============================   =================  =====================  =========  

Built-in Tools
--------------

Built-in tools are the scripts that can be run from command-line independently when you have installed the Cpipe package.


.. _Built-in tools:

============================   =====================  
Tool Name                      Modified from        
============================   =====================  
liftover
Venn Diagram
Conservation Plot
Correlation plot               bigwig_correlation
bamtofastq
BedClip
wigTobigwiggle
RegPotential
sample_contamination
============================   =====================  

.. _Bowtie:
.. _samtools:
.. _MACS2:
.. _MDSeqpos:
.. _BEDtools:
.. _External Tools:

External Tools
--------------

External Tools are the tools invoked by ChiLin by their path.

============================   =====================  ==================
Tool Name                      Download source         Version
============================   =====================  ==================
FastQC
R
Cython
MACS2                          `MACS site`_           2.0.10 20120605
CEAS                           `CEAS site`_           0.9.9.7
bedtools		       `bedtools site`_	      v2.16.2
pybedtools
samtools		       `SAMtools site`_	      0.1.17
Bowtie                         `Bowtie site`_         0.12.8
bedGraphToBigWig	       `UCSC utilities`_      v4
FastQC                         `FastQC site`_         v0.10.1
pdfTeX                         `pdfTex site`_         v1.40.10
IGV
============================   =====================  ==================    

