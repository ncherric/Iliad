..   Running Sphinx v5.0.2

.. _manual-main:

.. _installation: https://iliad.readthedocs.io/en/latest/getting_started/installation.html

===========================
ILIAD Genomic Data Pipeline
===========================

Iliad is a multi-functional **genomic data processing pipeline** implemented via the **Snakemake** workflow management system, Singularity or Docker container, 
and a handful of Conda environments.
Instances of any Singularity or Docker containers and Conda environments are automatically pulled down and 
created to build the right execution environment during workflow runtime, making light of the numerous required bioinformatic third-party software tools and dependencies.
The multi-functional feature is supported by a modularized Snakemake workflow design and includes independent Snakemake workflows to accommodate 
**many forms of genomic data** and translate them to comprehensible genotypes in a **variant call format (VCF) file**.

These are the currently supported genomic data starting points below:

* Raw Sequence Read Data (.fq)
* Stored Sequence Read Data (.cram)
* SNP Array Data (.idat)
* VCF data (.vcf/.vcf.gz)

Iliad is set apart by empowering genetic data management and processing for reference and target datasets independently or simultaneously.

What you need to do:

* **Follow the Installation and How-To Guides** for the appropriate module or submodule you are interested in running
* **Add your working directory PATH** to the top of the configuration file ``Iliad/config/config.yaml``
* **provide an FTP link** for the desired open-source raw or stored data into an Excel file found in the cloned ``Iliad/config`` directory OR **migrate locally hosted data** into the cloned Iliad repository

What Iliad will do:

* **Deliver the GWAS data** you need for your analyses as optimized VCF output

.. _main-getting-started:

---------------
Getting started
---------------

* Please read our paper for a full understanding of how Iliad can best serve you and your genetic data processing road blocks `replace with paper link here <https://github.com/ncherric/Iliad>`_.

* Quickly learn how to install Iliad - visit the Installation_ page

* Start processing your data right away by following one of the `HOW-TO GUIDES` on Iliad modules or `Submodules`

  * Modules

      *  Raw Sequence Read Data (.fq)

      *  Stored Sequence Read Data (.cram)

      *  SNP Array Data (.idat)


**Iliad Workflow Schematics**
****************************

.. image:: tutorial/img/Iliad-Figure1-forReadTheDocs.png
   :align: center
   :width: 600

.. image:: tutorial/img/Iliad-Figure2-forReadTheDocs.png
   :align: center
   :width: 600

.. toctree::
   :caption: Getting started
   :name: getting_started
   :hidden:
   :maxdepth: 1

   getting_started/installation
   getting_started/platform_preparation
   getting_started/demo
   getting_started/config

.. toctree::
   :caption: How-To Guides
   :name: tutorial
   :hidden:
   :maxdepth: 1

   tutorial/raw_sequence
   tutorial/stored_sequence
   tutorial/snp_array


.. toctree::
   :caption: Submodules
   :name: tutorial
   :hidden:
   :maxdepth: 1

   tutorial/lift_and_merge
   tutorial/liftover
   tutorial/merge_ref_target

.. toctree::
   :caption: Project Info
   :name: tutorial
   :hidden:
   :maxdepth: 1

   projectinfo/citing
   projectinfo/more_resources
   projectinfo/faq
   projectinfo/contributing
   projectinfo/credits
   projectinfo/changelog
   projectinfo/license

.. toctree::
   :maxdepth: 2
   :caption: Contents:
