..   Running Sphinx v5.0.2

.. _manual-main:

.. _installation: https://iliad-readthedocs.readthedocs.io/en/latest/getting_started/installation.html

===========================
ILIAD Genomic Data Pipeline
===========================

.. image:: https://img.shields.io/badge/snakemake-≥6.3.0-green.svg
    :target: https://snakemake.github.io

.. image:: https://img.shields.io/badge/python-≥3.8.0-brightgreen.svg
    :target: https://python.org

.. image:: https://img.shields.io/badge/Singularity-≥3.6.4-blue.svg
    :target: https://docs.sylabs.io/guides/3.6/user-guide/introduction.html

Iliad is a multi-functional **genomic data processing pipeline** implemented via the **Snakemake** workflow management system, singularity or docker container, and a handful of conda environments.
Instances of the Singularity or Docker container and Conda environments are automatically pulled down and created to build the right execution envioronment during workflow runtime, making light of the numerous required bioinformatic third-party software tools and dependencies.
The multi-functional feature is supported by a modularized Snakemake workflow design and includes independent snakemake workflows to accomodate **many forms of GENOMIC DATA** and translate them to comprehensible genotypes in a **variant call file (VCF)**.

* Raw Sequence Read Data (.fq)
* Stored Sequence Read Data (.cram)
* SNP Array Data (.idat)
* VCF data (.vcf/.vcf.gz)

Iliad is set apart by empowering genetic data management and processing for reference and target datasets, independently or simultaneously.

What you need to do:

* **provide an FTP link** for the desired open-source raw or stored data OR
* **migrate locally hosted data** into the cloned Iliad repository

What Iliad will do:

* deliver the 'GWAS data' you need for your analyses

.. _main-getting-started:

---------------
Getting started
---------------

* Please read our paper for a full understanding of how Iliad can best serve you and your genetic data processing road blocks `replace with paper link here <https://github.com/ncherric/Iliad-ReadtheDocs>`_.

* Quickly learn how to install Iliad - visit the Installation_ page

* Start processing your data right away by following along one of the `HOW-TO GUIDES` on Iliad modules

  *  Raw Sequence Read Data (.fq)
  *  Stored Sequence Read Data (.cram)
  *  SNP Array Data (.idat)

    .. make each of these into a link

.. toctree::
   :caption: Getting started
   :name: getting_started
   :hidden:
   :maxdepth: 1

   getting_started/installation

.. toctree::
   :caption: How-To Guides
   :name: tutorial
   :hidden:
   :maxdepth: 1

   tutorial/raw_sequence
   tutorial/stored_sequence
   tutorial/snp_array
   tutorial/merge_ref_target
   tutorial/liftover

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



.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
