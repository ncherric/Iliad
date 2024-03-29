
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org
.. _instructions: https://mamba.readthedocs.io/en/latest/installation.html
.. _platform: https://iliad-readthedocs.readthedocs.io/en/latest/getting_started/platform_preparation.html
.. _token: https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-https-urls
.. _creation: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-fine-grained-personal-access-token
.. _homepage: https://cloud.google.com/?hl=en

.. _getting_started/academic:

============================
Academic HPC cluster - Linux
============================

Step 1: Install Snakemake and Snakedeploy
*****************************************

Snakemake and Snakedeploy are best installed via the `Mamba package manager <https://github.com/mamba-org/mamba>`_ (a newer package manager for Conda replacement). 
From our testing, Conda is perfectly fine as an alternative, and you can use them interchangeably below.
If you have neither Conda nor Mamba, it can be installed via `Mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`_. 
For other options see `here <https://github.com/mamba-org/mamba>`_.

You will need to add Mamba to PATH so follow their instructions_ for such.

Here is my suggested install guidance for Conda on an academic HPC:

.. code-block:: console

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    $ CONDA_DIR=/opt/conda
    $ /bin/bash ~/miniconda.sh -b -p /opt/conda
    $ PATH=$CONDA_DIR/bin:$PATH
    $ conda init bash
    $ source ~/.bashrc

Given that Mamba (or Conda) is installed, run the following in your command line interface tool 
(i.e. `Putty <https://www.putty.org/>`_, 
`MacOS Terminal <https://support.apple.com/guide/terminal/open-or-quit-terminal-apd5265185d-f365-44cb-8b09-71a064a42125/mac>`_,
or `VS Code <https://code.visualstudio.com/>`_).
If you need guidance for a specific platform you are using, see the specific `platform`_ preparation guides and they will help lead you to this point.

.. code-block:: console

    $ conda create -c conda-forge -c bioconda --name iliadEnv snakemake=7.19.0 snakedeploy openpyxl pandas

to install Snakemake, Snakedeploy, and a necessary python library for parsing user input tables in an isolated environment.
For all following commands ensure that this environment is activated via


.. code-block:: console

    $ conda activate iliadEnv

**Estimated time: 15-20 minutes**

Step 2: Install Singularity
***************************

Typically users will not have root access on a shared HPC cluster in an academic setting. It is likely that Singularity is already installed on your school's system. 
In the case of Indiana University, a user can run the following command to load Singularity:

.. code-block:: console

    $ module load singularity/3.6.4

Your school will likely have a different version and you can check by either spamming the TAB button to autofill or with:

.. code-block:: console

    $ module avail

If your school does not have Singularity, you may have to ask your HPC IT team to install it or use an alternative Iliad installation method.

Step 3: Clone the Iliad repository and workflows
************************************************

.. Given that Snakemake and Snakedeploy are installed and available (see Step 1), the workflow can be deployed as follows.

First, create an appropriate project working directory (/path/to/project) on your system and enter it:

.. code-block:: console

    $ mkdir -p project
    $ cd project

In the next step, you will clone the Iliad repo. This will create an Iliad directory that you will cd into.
If you are not an active github user, you may have to create an account and a personal access token that is entered 
for password when prompted to do so on the command line in the following step. 
Here is a link for token_ information and creation_.

.. **OPTION 1: snakedeploy**

.. .. code-block:: console

..     $ snakedeploy https://github.com/ncherric/Iliad . --tag v1.0.0
..     $ cd Iliad

**Clone the repository using git:**

.. Git clone the `GitHub repository <https://github.com/ncherric/Iliad>`_.

.. code-block:: console

    $ git clone https://github.com/ncherric/Iliad.git
    $ cd Iliad

Now, there should be an ``./Iliad`` directory cloned into your ``/path/to/project/`` like such ``/path/to/project/Iliad/``.
And your current working directory should be ``/path/to/project/Iliad/``.

Two important folders found in the Iliad directory are **workflow** and **config**.
The ``workflow`` contains rules and scripts that a designated Snakefile in Iliad call on to run a specific module.
The ``config`` contains one configuration file ``Iliad/config/config.yaml`` which will be modified in the next step in order to configure the workflow to your needs.
It also contains ``Excel`` files and ``TSV`` files where you will input your sample information.

.. **side note**
.. ( Once this pipeline is publicly available, and added to the Snakemake Workflow Catalog, run below. For now, just **clone ABOVE** )


Step 4: Configure Workflow
**************************

There are 2 methods: Automatic and Manual 

**A) Automatic**

.. code-block:: console

    $ cd config
    $ python auto_config.py
    # Now you will see interactive prompts. If you want to follow default and tutorial, here are your options that you should copy and paste, individually when prompted.
    </PATH/TO/project/Iliad/>
    config/UserSampleTable.csv
    homo_sapiens
    104
    GRCh38
    $ mv modified.yaml config.yaml
    $ cd ..

There will be interactive questions on the command line that will ask you to enter your 1) working directory 2) sample table file with download URLs and 3) reference assembly info for download. 
Answer the interactive prompts accordingly and then press RETURN/ENTER.
NOTE: Using this command-line interactive prompt to update the config.yaml file will erase all comments and notes in your ``config.yaml`` file. 
There is an additional ``config-commented.yaml`` that you can refer to if you have questions about settings.


**B) Manual**

To configure this workflow, modify ``Iliad/config/config.yaml`` according to your needs. 
The file is clearly denoted into sections that you can change according to your needs. 
There are many defaults set that you do not have to change. The one most important change you will have to make is the following:

.. code-block:: console

    $ nano config/config.yaml

And INSERT your working directory path where NEED PATH HERE is. should look like this: **/path/to/project/Iliad/**

.. code-block:: yaml

    #####################################
    #####################################
    #####################################

    #  #  # USER INPUT VARIABLES  #  #  #

    #####################################
    #####################################
    #####################################

    # You must insert your /PATH/TO/Iliad/
    # use 'pwd' command to find your current working directory when you are inside of Iliad directory
    # e.g. /path/to/Iliad/ <---- must include forward slash at the end of working directory path

    # must include forward slash, '/', at the end of working directory path
    workdirPath: NEED PATH HERE


Step 5: Run workflow
********************

Given that the workflow has been properly deployed and configured, and your conda environment is activated, it can be executed as follows.

For running the workflow while deploying any necessary software via singularity and conda (using the Mamba package manager by default), run Snakemake with

.. code-block:: console

    $ snakemake -p --use-singularity --use-conda --cores 1 --jobs 1 --default-resource=mem_mb=10000 --latency-wait 120

Other workflows
********************

When executing the workflow, Snakemake will automatically find the main Snakefile in the workflow subfolder which is the **Raw Sequence Read Data** module.
However, there are other Snakefiles that are specific to the other Modules that you will have to call using ``--snakefile [desired module snakefile]``

* ``--snakefile workflow/Snakefile`` -> Raw Sequence Read Data
* ``--snakefile workflow/cram_Snakefile`` -> Stored Sequence Read Data
* ``--snakefile workflow/snpArray_Snakefile`` -> SNP Array Data
* ``--snakefile workflow/LiftoverTo38_Snakefile`` -> Submodule to liftover GRCh37 assembly VCF data to GRCh38 assembly
* ``--snakefile workflow/LiftoverTo37_Snakefile`` -> Submodule to liftover GRCh38 assembly VCF data to GRCh37 assembly
* ``--snakefile workflow/merger_Snakefile`` -> Submodule to merging list of VCFs
* ``--snakefile workflow/mergeRefTarget_Snakefile`` -> Submodule that will merge your processed Reference and Target data if you have previously completed both modules 

Visit the How-To Guides pages for further info about each of the Modules and Submodules in the bulleted list above.

This example bewlow is for the `Stored Sequence Read Data <https://iliad-readthedocs.readthedocs.io/en/latest/tutorial/stored_sequence.html>`_

.. code-block:: console

    $ snakemake -p --use-singularity --use-conda --cores 1 --jobs 1 --snakefile workflow/cram_Snakefile --default-resource=mem_mb=10000 --latency-wait 120

