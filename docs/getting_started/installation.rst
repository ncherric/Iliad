.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org
.. _instructions: https://mamba.readthedocs.io/en/latest/installation.html
.. _platform: https://iliad-readthedocs.readthedocs.io/en/latest/getting_started/platform_preparation.html
.. _token: https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-https-urls
.. _creation: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-fine-grained-personal-access-token
.. _login: https://cloud.google.com/?hl=en
.. _getting_started/installation:

============
Installation
============


Iliad is Snakemake workflow management system and can be cloned from the `GitHub repository <https://github.com/ncherric/Iliad>`_.
It largely reduces the learning curve and software dependency installations by taking advantage of pre-built execution environments for Iliad. 
It still requires the installation of the native package manager for the Snakemake application itself and Singularity.
Many HPC users at an academic institution will already have Singularity on their clusters.

.. image:: img/How-To-INTRO_Trim.mp4
   :align: center
   :width: 600

Google Cloud Platform server - Linux (Ubuntu)
====================================================================

First, you will have to use your google account or create one and login_.

Step 1: Install Snakemake and Snakedeploy
*****************************************

Snakemake and Snakedeploy are best installed via the `Mamba package manager <https://github.com/mamba-org/mamba>`_ (a newer package manager for Conda replacement). 
From our testing, Conda is perfectly fine as an alternative, and you can use them interchangeably below.
If you have neither Conda nor Mamba, it can be installed via `Mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`_. 
For other options see `here <https://github.com/mamba-org/mamba>`_.

You will need to add Mamba to PATH so follow their instructions_ for such.

Here is my suggested install guidance for Conda:

.. code-block:: console

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    $ CONDA_DIR=/opt/conda
    $ sudo /bin/bash ~/miniconda.sh -b -p /opt/conda
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

Follow the `instructions guide <https://docs.sylabs.io/guides/3.6/user-guide/quick_start.html>`_ provided by Singularity to install the program under the 
Operating System of your choosing.
If you are on an HPC server, it may have a Singularity module available to load already.

Although it is highly recommended to follow the Singularity specific installation, we have some install advice that worked for us and may work best for you.
Again, please see our specific `platform`_ preparation guides to set up your machine to this point if you are on a Windows or Mac machine. 
Singularity cooperates best with a Linux OS. It requires GO, so the following code blocks incorporate this dependency.

You might need to install system dependencies.. For example, if you just created a VM Instance on Google Cloud Platform:

Debian-based systems, such as Ubuntu 20.04 (x86_64)

Ensure repositories are up-to-date and install debian packages for dependencies

.. code-block:: console

	$ sudo apt-get update

	$ sudo apt-get install -y \
	build-essential \
	libseccomp-dev \
	pkg-config \
	squashfs-tools \
	cryptsetup \
	curl wget git

CentOS/RHEL systems:

Install basic tools for compiling, ensure EPEL repository is available, install RPM packages for dependencies

.. code-block:: console

	$ sudo yum groupinstall -y 'Development Tools'
    
	$ sudo yum install -y epel-release

	$ sudo yum install -y \
	libseccomp-devel \
	squashfs-tools \
	cryptsetup \
	wget git

Install GO and put it in your PATH

.. code-block:: console

	$ export GOVERSION=1.18.1 OS=linux ARCH=amd64  # change these variables as you need
	$ wget -O /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz https://dl.google.com/go/go${GOVERSION}.${OS}-${ARCH}.tar.gz
	$ sudo tar -C /usr/local -xzf /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz
	$ echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
	$ source ~/.bashrc

Install Singularity

.. code-block:: console

	$ export VERSION=3.8.7  # adjust this as necessary
	$ wget https://github.com/apptainer/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
	$ tar -xzf singularity-${VERSION}.tar.gz
	$ cd singularity-${VERSION}

.. code-block:: console

	$ ./mconfig
	$ make -C ./builddir
	$ sudo make -C ./builddir install

Make sure Singularity is in your PATH and accessible by Snakemake.

.. code-block:: console

	$ singularity version

Edit the Singularity configuration file to allow SHARED LOOP DEVICES.

.. code-block:: console
   
	$ sudo nano /usr/local/etc/singularity/singularity.conf

Change the SHARED LOOP DEVICES to yes. "shared loop devices = yes"

.. code-block:: yaml

    # SHARED LOOP DEVICES: [BOOL]
    # DEFAULT: no
    # Allow to share same images associated with loop devices to minimize loop
    # usage and optimize kernel cache (useful for MPI)
    shared loop devices = yes

Return to Home directory

.. code-block:: console

	$ cd ~

**Estimated time: 15-60 minutes**

Step 3: Clone the Iliad repository and workflows
************************************************

.. Given that Snakemake and Snakedeploy are installed and available (see Step 1), the workflow can be deployed as follows.

First, create an appropriate project working directory (/path/to/project-workdir) on your system and enter it:

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

Now, there should be an ``./Iliad`` directory cloned into your ``/path/to/project-workdir/`` like such ``/path/to/project-workdir/Iliad/``.
And your current working directory should be ``/path/to/project-workdir/Iliad/``.

Two important folders found in the Iliad directory are **workflow** and **config**.
The ``workflow`` contains rules and scripts that a designated Snakefile in Iliad call on to run a specific module.
The ``config`` contains one configuration file ``Iliad/config/config.yaml`` which will be modified in the next step in order to configure the workflow to your needs.
It also contains ``Excel`` files and ``TSV`` files where you will input your sample information.

.. **side note**
.. ( Once this pipeline is publicly available, and added to the Snakemake Workflow Catalog, run below. For now, just **clone ABOVE** )




Step 4: Configure Workflow
**************************

To configure this workflow, modify ``Iliad/config/config.yaml`` according to your needs. 
The file is clearly denoted into sections that you can change according to your needs. 
There are many defaults set that you do not have to change. The one most important change you will have to make is the following:

.. code-block:: console

    $ nano config/config.yaml

.. code-block:: console

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
    workdirPath: /Insert/path/to/Iliad/


Step 5: Run workflow
********************

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via singularity and conda (using the Mamba package manager by default), run Snakemake with

.. code-block:: console

    $ snakemake -p --use-singularity --use-conda --cores 1 --jobs 1 --default-resource=mem_mb=10000 --latency-wait 120


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
