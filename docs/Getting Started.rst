===============
Getting Started
===============


Setting up ChiLin
===================


install quick-start
----------------------
In order to make installation easier, we write a script, download.sh, run with help instruction


test installation
---------------------
we sample down 1000 reads for testing the ChiLin installation state
the test files are in the *tests* folder

 Config file

 To run the tests, you need simply need to run::

    Chilin.py gen -s hg19

 then, replace the *ChiLinjinja.conf basis* section with test.conf
    *basis* section, run::
    ChiLin.py run -c ChiLinjinja.conf -t TF


Synopsis
--------------

.. envvar:: [meta]

    Lists all the meta-data of current workflow.
    
    Consist of the following options:

    .. envvar:: dataset.ID

        The name for the dataset, which will be the value of :envvar:`${DatasetID}`
	
	Limit: a string (1) consist of ``numbers``, ``alphabets`` or ``'_'`` (2) shorter than 20 characters

    .. envvar:: species
        The name of species, written to the QCreport and log

	Limit: a string (1) consist of ``numbers``, ``alphabets`` or ``'_'`` (2) shorter than 20 characters

    .. envvar:: factor
        
        The name of species, writen to DC summary and QCreport, log 
	Limit: a string (1) come from GO standard term
     
    .. envvar:: TF or Histone

        The judgement of the factor for MACS shift-size choice and assess motif QC measurement 

    Limit: Logic value, True for TF and False for Histone

    .. envvar:: assembly

        The assembly version, written to the QCreport and log

	Limit: a string (1) consist of ``numbers``, ``alphabets`` (2) shorter than 10 characters

    .. envvar:: treatment

       The paths of treatment files

       Limit: absolute or relative ``path`` of files in :ref:`supported formats<raw data>`

    .. envvar:: control

       The paths of treatment files

       Limit: absolute or relative ``path`` of files in :ref:`supported formats<raw data>`
       
.. envvar:: [ext]

    The external data and external tools to use. Read :ref:`External Data` and  :ref:`External Tools` for a full explanation.
    
.. envvar:: [steps]

    meta

See :ref:`simpest_config` for a quick view of how to use this pipeline.

.. _simpest_config:

Simpest config
-----------------

Here is one of the simpest Cpipe workflow you can make.

.. literalinclude:: demo/start_chilin.conf
      :language: ini
      :linenos:

Use your own path of :ref:`Raw Data<Raw Data>` to replace the Line 6. And use the path of the directory used to store :ref:`External Data` to replace Line 9.

When saved to ``starter_chilin.conf``, this config file can construct a powerful pipeline via:

::

    $ chilin hello_cpipe.conf

When it finished about 2 hours later, you will get :ref:`Processed Data<Processed Data>` and a :ref:`PDF report`.


    
Examples about replicates
=========================

For example, there are five raw files. Three of them are replicates for ``treatment`` and two of two for ``control``.

The file names may look like:

::

    demo_treat1.fastq
    demo_treat2.fastq
    demo_treat3.fastq
    demo_control1.fastq
    demo_control2.fastq


Then you can write the :envvar:`[meta]` section like this:

.. code-block:: ini
   :linenos:
    [UserInfo]
    User = testuser
    datasetID = testid
    species = hg19
    factor = testfactor
    treatpath = /mnt/Storage/home/qinq/treat1.fastq,/mnt/Storage/home/qinq/treat2.fastq
    controlpath = /mnt/Storage/home/qinq/control1test.fastq
    OutputDirectory = /mnt/Storage/home/qinq/testchilin3
    ...

Replace the commented in Line 2, Line 3 and Line 4 and complete other sections. Then load it with Cpipe.

For the notation of output files, the :envvar:`${DatasetID}` will be ``demo_replicate``. The :envvar:`${treat_rep}` will be ``1``, ``2`` and ``3``. The :envvar:`${control_rep}` will be ``1`` and ``2``.

Common usage of ChiLin modules
===============================
Analysis Pipeline
-------------------

Quality Control Pipeline
-------------------------

Analyzing single datasets
---------------------------


Analyzing datasets in batch
-------------------------------





*Getting help*
------------------
Any question on installation or runnning is appreciated, please mail
to qinqianhappy@gmail.com.


