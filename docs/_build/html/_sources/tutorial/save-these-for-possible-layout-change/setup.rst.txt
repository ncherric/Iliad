.. _snpArrayTutorial-Basics:

.. _MEGA: https://support.illumina.com/array/array_kits/infinium-multi-ethnic-global-8-kit.html
.. _installation: https://iliad-readthedocs.readthedocs.io/en/latest/getting_started/installation.html
.. _module: raw sequence read module here

Setup
::::::::::

Once the Installation_ of Iliad and its two dependencies has been completed, 
you will find your new working directory within the ``PATH/TO/Iliad`` folder.
Make sure your current working directory is in this cloned repo as stated in the installation.

.. code:: console

   $ cd Iliad

In that working directory you will find there are a number of directories with files and code to run each of the module pipelines.

**FIRST**, there is a ``data/snp_array/idat`` directory with a ``readme.md`` file. You must place all of your ``.idat`` files in this folder.
There should be two ``.idat`` files for each sample: one green ``_Grn.idat`` and one red ``_Red.idat``. 

**SECOND**, there is a configuration file with some default parameters, however, you MUST at least change the ``workdirPath`` parameter to the appropriate 
path leading up to and including ``/Iliad`` e.g. ``/Path/To/Iliad/``. The configuration file is found in ``config/config.yaml``.

.. code:: python

    workdirPath: /Path/To/Iliad/

Some other parameters that are pre-set and you might consider changing to your project needs include:

* Homo sapiens GRCh38 release 104 reference genome

.. code:: python

    ref:
      species: homo_sapiens
      release: 104
      build: GRCh38

* Illumina MEGA_ microarray GRCh38 support and product files

.. code:: python

    urlProductFiles:
      manifest:  https://webdata.illumina.com/downloads/productfiles/multiethnic-global-8/v1-0/build38/multi-ethnic-global-8-d2-bpm.zip
      mzip: multi-ethnic-global-8-d2-bpm.zip
      cluster: https://webdata.illumina.com/downloads/productfiles/multiethnic-global-8/v1-0/infinium-multi-ethnic-global-8-d1-cluster-file.zip
      czip: infinium-multi-ethnic-global-8-d1-cluster-file.zip
    urlSupportFiles:
      physicalGeneticCoordinates: https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/multiethnic-global/multi-ethnic-global-8-d2-physical-genetic-coordinates.zip
      pzip: multi-ethnic-global-8-d2-physical-genetic-coordinates.zip
      rsidConversion: https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/multiethnic-global/multi-ethnic-global-8-d2-b150-rsids.zip
      rzip: multi-ethnic-global-8-d2-b150-rsids.zip
      rfile: Multi-EthnicGlobal_D2_b150_rsids.txt

**THIRD**,
each module pipeline has a specific ``Snakefile``.
Snakemake will automatically detect the main snakefile, which is named excatly as such and found in the ``workflow`` directory: ``workflow/Snakefile``.
Iliad reserves the main snakefile for the main module, specifically the raw sequence read data module_.
This means the user must specify which ``Snakefile`` will be invoked with 

.. code:: console

    $ snakemake --snakefile workflow/snpArray_Snakefile

and combined with other user-specified snakemake flags, of course, like ``--cores``.

In this module, the SNP Array Snakefile is also located in the ``workflow`` directory: ``workflow/snpArray_Snakefile``.
Users must invoke this ``snpArray_Snakefile`` in order to run the workflow for This **SNP ARRAY MODULE**.




