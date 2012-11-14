ChiLin
------


ChiLin is a pipeline for ChIP-seq data written in pure Python. It provides workflow of mapping(Bowtie),
peak calling(MACS2) and motif finding(MDSeqpos). Moreover, it provides quality
control tools for reads, replication and peaks. It's intended to do more things
with less configuration.

Nutshell
--------

Here is a workflow of using ChiLin::


     1. type in this command to generate a config file:

     ChiLin gen -s hg19 > example.conf

     2.


     [meta]
     user =
     dataset_id =
     species = hg19
     factor =
     treatments =
     controls =
     output_dir =


