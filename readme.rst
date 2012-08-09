
About ChiLin
============

ChiLin is the Chinese Unicorn. It is sometimes called the dragon horse, a mythical creature that has a lot of legends and stories associated to it.

We use ChiLin for the name of our project as it has similar alphabet combination with **ChI**\ P-seq pipe\ **lin**\ e.

Usage of ChiLin
===============

usage: ChiLin.py [-h] {run,gen} ...

ChiLin : A clear ChIP-seq pipeline

positional arguments:
  {run,gen}   sub-command help
    gen       generate a template of config file
    run       run pipeline using a config file

optional arguments:
  -h, --help  show this help message and exit



ChiLin-run: Run ChiLin pipeline using a config file
---------------------------------------------------

usage: ChiLin.py run [-h] -c CONFIG -t {Dnase,Histone,TF} [-m {mean}] [-p TOP_PEAKS] [--threads {1,2,3,4,5,6,7,8}] [-d STEP_END] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  
  -c CONFIG, --config CONFIG
                        specify the config file to use
			
  -t {Dnase,Histone,TF}
                        the most important option for ChiLin specify the analysis type and the shiftsize {Dnase: 50, Histone and TF:73} of MACS2
			
  -m {mean}             specify method for correlation plot
  
  -p TOP_PEAKS          specify peaks number for CEAS
  
  --threads {1,2,3,4,5,6,7,8}
                        How many threads can be used
			
  -d STEP_END           specify the end step of pipeline, 1 for bowtie, 2 for macs, 3 for venn and correlation, 4 for ceas, 5 for conservation, 6 for motif, Note: if you only have bed file, start from 2
  
  --debug               debug mode




ChiLin-gen: A config template generator for ChiLin
--------------------------------------------------

usage: ChiLin.py gen [-h] --species {hg19,mm9}

optional arguments:
  -h, --help            show this help message and exit
  
  --species {hg19,mm9}  

Quick-install
=============

R
-

install R, recommended the newest version
in R session
> install.packages('gplots')
> install.packages('Rcolorbrew')


Python
------

In the terminal, change directory to the folder 
python setup.py install # need root

External program
----------------

We created a package containing all the sources
