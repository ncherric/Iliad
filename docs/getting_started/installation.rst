.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org
.. _instructions: https://mamba.readthedocs.io/en/latest/installation.html
.. _platform: https://iliad-readthedocs.readthedocs.io/en/latest/getting_started/platform_preparation.html
.. _token: https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-https-urls
.. _creation: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-fine-grained-personal-access-token
.. _homepage: https://cloud.google.com/?hl=en
.. _Docker: https://hub.docker.com/repository/docker/ncherric/iliad/general
.. _Singularity: https://cloud.sylabs.io/library/ncherric/iliad/igdp-container

.. _getting_started/installation:

==============================
Installation and Configuration
==============================


Iliad is Snakemake workflow management system and can be cloned from the `GitHub repository <https://github.com/ncherric/Iliad>`_.
It largely reduces the learning curve and software dependency installations by taking advantage of pre-built execution environments for Iliad. 
It still requires the installation of the native package manager for the Snakemake application itself and Singularity with certain machine configurations.
Many HPC users at an academic institution will already have Singularity on their clusters just be sure to load the module for Snakemake to find the Iliad Singularity_ Image. 
The Iliad Docker_ Image provides an alternative method to using Singularity. 

https://youtu.be/9CCnaLlUFG4

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/Wu0EdBP_CD0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

Quick intro about Iliad

Google Cloud Platform server - Linux
====================================================================

LINK

Academic HPC cluster - Linux
============================

LINK


Local machine - Docker
====================================================================

LINK