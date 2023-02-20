.. _tutorial-snpArrayTutorial:

========================
Iliad SNP Array Tutorial
========================

.. hyperlinks
.. _Iliad: https://iliad-readthedocs.readthedocs.io/
.. _Snakemake: https://snakemake.readthedocs.io
.. _Illumina: https://support.illumina.com/
.. _MEGA: https://support.illumina.com/array/array_kits/infinium-multi-ethnic-global-8-kit.html
.. _GCP: https://cloud.google.com/
.. _iaap-cli: https://support.illumina.com/downloads/iaap-genotyping-cli.html
.. _EULA: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/iaap/Illumina%20Array%20Analysis%20Platform%20IAAP%201.1%20EULA.pdf
.. _gtc2vcf: https://github.com/freeseek/gtc2vcf
.. _bcftools: https://samtools.github.io/bcftools/bcftools.html
.. _installation: https://iliad-readthedocs.readthedocs.io/en/latest/getting_started/installation.html
.. _module: raw sequence read module here
.. _slides: https://slides.com/johanneskoester/snakemake-tutorial


This tutorial introduces the genome-wide SNP array data processing module of the Iliad_ workflow developed using Snakemake workflow language.
Please visit Snakemake_ for specific details. In general, though, each module is composed of rules. These rules define how output files are generated from input files while 
automatically determining dependencies amongst the rules. A ``DAG`` (directed acyclic graph) of jobs will be built each time to account for all of the samples and jobs 
that will executed either via job scheduler or local cores and will execute in parallel if multiple jobs are declared.
Because of the Snakemake workflow system design, the **Iliad** workflow is scalable from single core machines to HPC clusters with job schedulers.

The **SNP array module** is designed to process target data in your lab. Iliad is currently limited to Illumina_ microarray raw data processing and is configured for the 
human genotyping Infinium Multi-Ethic Global-8 Kit (MEGA_).
We ensured no bioinformatics knowledge is needed to run this module with the help of external test runs performed on Google Cloud Platform (GCP_).


**SNP Array Module DAG**
.. image:: img/snp_array_module_dag.png
   :align: center

.. toctree::
   :maxdepth: 2

   background
   basics
   setup
   
