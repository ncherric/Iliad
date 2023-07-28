<snippet>
  <content>

<p align="center">
  <img width="496"  src="docs/readme-img/ILIAD-logo.png" alt="Iliad logo">
</p>


[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-green.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-≥3.8.0-brightgreen.svg)](https://python.org)
[![Singularity](https://img.shields.io/badge/Singularity-≥3.6.4-blue.svg)](https://docs.sylabs.io/guides/3.6/user-guide/introduction.html)
---------

# Iliad's workflows

## Main Modules
<p align="center">
  <img width="496"  src="https://github.com/ncherric/Iliad/blob/main/docs/readme-img/Iliad-Figure1-forReadTheDocs.png?raw=true">
</p>

## A comprehensive Submodule
<p align="center">
  <img width="496"  src="https://github.com/ncherric/Iliad/blob/main/docs/readme-img/Iliad-Figure2-forReadTheDocs.png?raw=true">
</p>


Iliad is a suite of Snakemake workflows for genomic data processing. The modularization feature of Snakemake accomodates raw sequence data (.fq), stored sequence data (.cram), Illumina microarray data (.idat), and variant call format data (.vcf).

---------
## Usage

Please visit our documentation hosted on Read The Docs: https://iliad.readthedocs.io/

## Intro Page

ILIAD 

Iliad is a multi-functional toolkit comprised of **genomic data processing pipelines** implemented via the **Snakemake** workflow management system, Singularity or Docker container, 
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

## Getting started

* Please read our paper for a full understanding of how Iliad can best serve you and your genetic data processing road blocks: 

* Quickly learn how to install Iliad - visit the Installation page at  https://iliad.readthedocs.io/

* Start processing your data right away by following one of the `HOW-TO GUIDES` on Iliad modules or Submodules

