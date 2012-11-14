=========================
INSTALL Guide for ChiLin
=========================


R packages
----------

install R, recommended the newest version
in R session
> install.packages('gplots')
> install.packages('RColorBrewer')

External tools
--------------
Before you run, please modify ChiLinjinja.conf.sample to up to your dependency.
or run download.sh to automatically config as ChiLinjinja.conf.new
We created a bash shell script `download.sh` to prepare all the dependent data and program for you.
# cd path_to_chilin
# sudo bash install.sh -p download_data_path # should mkdir under the chilin directory

ChiLin
------
In the terminal, change directory to the folder
#python setup.py install # need root