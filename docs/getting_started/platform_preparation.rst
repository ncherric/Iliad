.. _Installation: https://mamba.readthedocs.io/en/latest/installation.html

.. _getting_started/platform_preparation:

====================
Platform Preparation
====================
########################

Linux
======

If you are on a High-Performance Cluster (HPC) or remote server that runs on a Linux OS, like Ubuntu or Red Hat for example, 
then please follow the instructions on the previous Installation_ page!

|
|

Windows / PC
============

* Step 1. Freely download the UBUNTU app from the Windows store

* Step 2. Open the UBUNTU app and create a username and pw

* Step 3. Run the following:

obtain mamba package manager:

.. code-block:: console

    $ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"

install mambaforge

.. code-block:: console

    $ bash Mambaforge-$(uname)-$(uname -m).sh -b

* Step 4. RESTART

* Step 5. Install Singularity following the instructions on the Installation_ page.

* Step 6. Install conda environment for Iliad which contains Snakemake

.. code-block:: console

    $ conda create -c conda-forge -c bioconda --name iliadEnv snakemake=7.19.0 snakedeploy openpyxl pandas

* Step 7. Change directories in the UBUNTU command line to get into proper location of your project folder

.. code-block:: console

    $ cd mnt/DRIVE/Path/To/Iliad/

* Step 8. Activate your conda environment (IliadEnv)

.. code-block:: console

    $ Mamba activate IliadEnv

* Step 9. Get your current location to add to the configuration file

.. code-block:: console

    $ pwd

* Step 10. Copy that location and paste into the ``Iliad/config/config.yaml`` as your working directory. This should start with ``/mnt/`` and end with a forward slash ``/`` like so  ``/mnt/Drive/Path/To/Iliad/``

|
|

MacOS
=====

Docker container
****************


Oracle VM
*********

* Step 1. Download the latest version of `Oracle VM VirtualBox <"https://www.oracle.com/virtualization/technologies/vm/downloads/virtualbox-downloads.html#vbox">_`. 

* Step 2. Use the GUI - there is an easy install that makes editing a virtual machine in the mac quite simple. This required since Singularity needs to run on Linux.

* Step 3. Use a pre-made singularity container on sylabs

.. code-block:: console

    export VM=sylabs/singularity-ce-3.8-ubuntu-bionic64 && \
    vagrant init $VM && \
    vagrant up && \
    vagrant ssh

* Step 4. Gain control through the virtal machine application. This removes the need to go editing the vagrantfile.

* Step 5. Increase your storage space alotted for the virtual machine.
