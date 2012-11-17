==============================
Advanced sections for ChiLin
==============================
This sections for developers who are interested in our development
and want to participate
**Developers**: in order to be sure they are developing the right project that fulfills requirements provided in this document.
**Users**: in order to get familiar with the idea of the project and suggest other features that would make it even more functional. 

Dive into options
===================

step control
---------------
If you don't want to go from beginning to the end of the pipeline, say, you just want the peaks calling results,
the choice is the *-d* option, you could end at wherever you want.

Debug mode
-------------
using *--debug*
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

if further analysis requirements is needed, you may want to ajust the default options we set here.
Dive into conf
===============

for ceas
-------------------
* bipromotor sizes
* promotor sizes

* rel_dist ::


for conservation
------------------
* w width of conservation plot region


for macs2
-----------------
* --shift-size --nomodel: optional



In-depth Configuration
------------------------

Dive into rule
===============
We separate Name Rules for output and temporary files from analysis codes part for easier to maintain,
If you don't feel comfortable of our name ways, it's simple for you to adjust it.

