ChiLin
------


ChiLin is a pipeline for ChIP-seq data written in pure Python. It provides workflow of mapping(Bowtie),
peak calling(MACS2), Cis-elementary visualization(CEAS) and motif finding(MDSeqpos). Moreover, it provides quality
control tools for reads, replication and peaks. It's intended to do more things
with less configuration.

Note
----

This is an early-stage snapshot of ChiLin. To get the latest version of Chilin, please go to [this repository](https://github.com/cfce/chilin)


Nutshell
--------

Here is a workflow of using ChiLin::


     1. type in this command to generate a config file:

     ChiLin gen -s hg19 > a_meaningful_name.conf

     2. modify the basis part to meet your requirement
        [basis]
        user = foo
        id = 2012
        time = 2102xxx
        species = {{ species }}
        factor = ESR1
        treat = absolutepath
        control = absolutepath
        output = absolutepath
     3. run the analysis part
        ChiLin run -c a_meaningful_name.conf -t TF(orHistone)
