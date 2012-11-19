===============
Getting Started
===============

Setting up ChiLin
===================

install quick-start
----------------------
In order to make installation easier, we write a script,
`download.sh`, run with help instruction. This script could get access
to the required data, such as bowtie index for hg19 and mm9, and all
stable dependent program.
basicly using the command below::

test installation
---------------------
we sample down 1000 reads for testing the ChiLin installation state
the test files are in the *tests* folder

- Generating Configuration file
 To run the tests, you need simply need to run::

    Chilin.py gen -s hg19

.. literalinclude:: ../chilin/tests/test.conf
      :language: ini
      :linenos:

 then, replace the *ChiLin.conf basis* section with
 chilin/tests/test.conf *basis* section. Fill in the treat and control in absolute
 path in the chilin/tests/ directory, the data are sampled down for
 speed test.

- run::

    ChiLin.py run -c ChiLin.conf -t TF

This should be run all dependent program with dependent data.

Common usage of ChiLin modules
===============================

Analyzing single datasets
---------------------------
Here we recommend a dataset from `GSE30624` for testing whole pipeline
with :ref:`QC report` generated.



Simpest config
--------------------------

Here is one of the simpest ChiLin_ workflow you can make.

.. literalinclude:: ../chilin/tests/test.conf
      :language: ini
      :linenos:

Use your own path of :ref:`Raw Data<Raw Data>` to replace the Line 6. And use the path of the directory used to store :ref:`External Data` to replace Line 9.
When saved to ``test.conf``, this config file can construct a powerful pipeline via:

::

    > ChiLin run -c test.conf -t TF --debug

When it finished about a few minutes later, you will get :ref:`Processed Data<Processed Data>` and a :ref:`PDF report`.
Replace the commented in Line 2, Line 3 and Line 4 and complete other sections. Then load it with Cpipe.
For the notation of output files, the :envvar:`${DatasetID}` will be ``demo_replicate``. The :envvar:`${treat_rep}` will be ``1``, ``2`` and ``3``. The :envvar:`${control_rep}` will be ``1`` and ``2``.


Analyzing datasets in batch
-------------------------------
For *autodc* module for running ChiLin, we provide a fork version
chilin-mvc_ with tools for auto downloading GEO data and analyzing on
server.

*Getting help*
==============
Any question on installation or runnning is appreciated, please mail
to qinqianhappy@gmail.com.


.. chilin-mvc:: https://bitbucket.org/Alvin_Qin/chilin-mvc
.. ChiLin:: https://bitbucket.org/shenglinmei/chilin
