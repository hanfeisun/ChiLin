
This page is intended for:

**Developers**: in order to make sure they're using the right format of data and right version of tool to test

**Users**: in order to know where they should go to download these data and tools


====
Data
====

Built-in Data
-------------

The Cpipe package includes all the build-in data for hg19 and mm9. For other species, you may need to download these data from data source or custom it yourself.
 
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

=====
Tools
=====

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


External Tools are the tools invoked by Cpipe by their path.

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


.. _MACS site: https://github.com/taoliu/MACS
.. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. _MDSeqPos site: https://bitbucket.org/cistrome/cistrome-applications-harvard/src/c477732c5c88/mdseqpos
.. _bedtools site: http://code.google.com/p/bedtools/
.. _SAMtools site: http://samtools.sourceforge.net/
.. _Bowtie site: http://bowtie-bio.sourceforge.net/index.shtml
.. _UCSC utilities: http://hgdownload.cse.ucsc.edu/admin/exe/
.. _UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables
.. _Cistrome site: http://cistrome.org/~hanfei
.. _FastQC site: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _pdfTex site: http://www.tug.org/applications/pdftex/ 

========
Workflow
========

.. digraph:: main
    
    rankdir=TB
    size="20,15"


.. digraph:: foo

    rankdir=TB
    size="15,15"
    edge[arrowhead=open]

    start[shape=circle, label="", style=filled]
    end[shape=doublecircle, label="", style=filled]

    readconf[shape=box,style=rounded, label="Read config"]
    bowtie[shape=box,style=rounded, label="Run Bowtie"]
    rawQC[shape=box,style=rounded, label="Run RawQC"]
    mappingQC[shape=box,style=rounded, label="Run MappingQC"]
    macs2[shape=box,style=rounded, label="Run MACS2"]
    peakcallingQC[shape=box,style=rounded, label="Run PeakcallingQC"]
    ceas_seqpos[shape=box,style=rounded, label="Run CEAS/Seqpos"]
    venn[shape=box,style=rounded, label="Draw VennDiagram"]
    conservation[shape=box,style=rounded, label="Draw ConservationPlot"]
    annotationQC[shape=box,style=rounded, label="Run AnnotationQC"]

    
    ifmapped[shape=diamond, label="Mapped?"]
    ifrep[shape=diamond, label="Replicate?"]
    
    start -> readconf
    readconf -> rawQC
    rawQC -> ifmapped[headport=n, color="grey"]
    ifmapped -> mappingQC[label="[Yes]" tailport=s]
    ifmapped -> bowtie[taillabel="[No]" tailport=e]
    bowtie -> mappingQC
    mappingQC -> macs2[color="grey"]
    macs2 -> ifrep
    peakcallingQC -> ceas_seqpos[color="grey"]
    ifrep -> venn[label="[Yes]" tailport=s]
    ifrep -> conservation[label="[No]" tailport=e]
    venn -> conservation
    conservation -> peakcallingQC
    ceas_seqpos -> annotationQC
    annotationQC -> end[taillabel="Output Report"]

