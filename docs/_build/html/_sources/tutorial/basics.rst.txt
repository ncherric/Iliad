.. _snpArrayTutorial-Basics:

.. _iaap-cli: https://support.illumina.com/downloads/iaap-genotyping-cli.html
.. _EULA: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/iaap/Illumina%20Array%20Analysis%20Platform%20IAAP%201.1%20EULA.pdf
.. _gtc2vcf: https://github.com/freeseek/gtc2vcf
.. _bcftools: https://samtools.github.io/bcftools/bcftools.html

Basics
::::::::::

The raw files from an Illumina sequencer are bead array files found in raw intensity data ``.idat`` format.
These ``.idat`` files are to be converted into Genotype Call ``.gtc`` files using iaap-cli_ software. This software does have an 
end-user license agreement (EULA_) and is not included or distributed by Iliad. If the user chooses to configure a download of the 
program, it will be downloaded, independent from the Iliad repository distribution. 

The ``.gtc`` files are converted to a ``.vcf`` using bcftools_ plug-in gtc2vcf_. 
This requires a reference genome assembly and Iliad downloads the user-configured reference genome fasta files. 
Iliad is configured to download *Homo sapiens* GRCh38 release 104 as default.

Processing the ``.vcf`` is critical for realistic genetic data compatibility.
Custom python scripts are called in rules to rename unconventional loci names to standardized ``rs IDs`` using dbSNP files.
The default configuration file is set to download human_9606_b151_GRCh38p7 ``All_20180418.vcf.gz``.

Once the ``vcf`` is processed, quality controls are performed based on the GenTrain and ClusterSep scores.
Default thresholds, along with other default workflow configurations, can be found in your path to the configuration file: ``config/config.yaml``.
