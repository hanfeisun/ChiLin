
Workflow
==========

.. digraph:: foo

    rankdir=TB
    size="15,15"
    edge[fontsize="11" arrowhead=open]

    start[shape=circle, label="", style=filled]
    end[shape=doublecircle, label="", style=filled]

    readconf[shape=box,style=rounded, label="Read config"]
    bowtie[shape=box,style=rounded, label="Run Bowtie"]
    rawQC[shape=box,style=rounded, label="Run RawQC"]
    mappingQC[shape=box,style=rounded, label="Run MappingQC"]
    macs2[shape=box,style=rounded, label="Run MACS2"]
    peakcallingQC[shape=box,style=rounded, label="Run PeakcallingQC"]
    ceas[shape=box,style=rounded, label="Run CEAS"]
    venn[shape=box,style=rounded, label="Draw VennDiagram"]
    conservation[shape=box,style=rounded, label="Draw ConservationPlot"]
    seqpos[shape=box,style=rounded, label="SeqPos"]
    motifQC[shape=box,style=rounded, label="MotifQC"]
    
    ifmapped[shape=diamond, label="Mapped?"]
    ifrep[shape=diamond, label="Replicate?"]
    
    start -> readconf
    readconf -> rawQC
    rawQC -> ifmapped[headport=n]
    ifmapped -> mappingQC[label="[Yes]" tailport=s]
    ifmapped -> bowtie[taillabel="[No]" tailport=e]
    bowtie -> mappingQC
    mappingQC -> macs2
    macs2 -> peakcallingQC
    peakcallingQC -> ceas
    ceas -> ceasQC
    ceasQC -> ifrep
    ifrep -> venn[label="[Yes]" tailport=s]
    ifrep -> conservation[label="[No]" tailport=e]
    venn -> conservation
    conservation -> conservQC
    conservQC -> seqpos
    seqpos -> motifQC
    
    motifQC -> end[taillabel="Output Report"]
    

