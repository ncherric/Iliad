.. _projectinfo/comingsoon:

===========
Coming Soon
===========


Sequence Read Archive retrieval
===============================

This feature will give users an option to automatically perform data downloads using the SRAtoolkit and access data from the NCBI Sequence Read Archive (SRA).


CRAM remapping
==============

Genomic data is commonly stored in CRAM format. Of course, these data files have already been mapped to a reference genome. 
To our knowledge, there is not a direct option for remapping a CRAM file, however, reverting a CRAM sample to a FASTQ or FASTA file is possible via samtools. 
This feature will include that step and steps thereafter for re-alignment to a reference genome of the user's choice. 


Variant discovery
=================
    
The current format of Iliad Sequence processing entails the use of variant genotyping using annotation files from the NYGC. We understand the importance of variant discovery 
and aim to streamline this option in a future release. Currently, you can manually (and quite easily) set this up by commenting out three lines of code and adding one new line
of code to ``./Iliad/workflow/Snakefile`` (for raw sequence reads) or ``./Iliad/workflow/cram_Snakefile`` (for cram sequence data) while we develop this method further. 
See the code blocks and images below to manually perform this.

**Raw sequence reads**

.. image:: img/coming-soon-Seq.png
   :align: center
   :width: 800

Raw sequence reads - file found in ``./Iliad/workflow/Snakefile``

You can delete the contents of ``./Iliad/workflow/Snakefile`` and copy/paste below to replace it:
.. code-block:: console

    singularity: "library://ncherric/iliad/igdp-container:v1.16"

    include: "rules/common.smk"
    #splits=config['NYGC']['numberOfSplitRegionsFiles']

    ##### Target rules #####
    rule all:
        input: 
            expand("results/vcf/Merged-chr{chroms}-allSamples-with-rsIDs.vcf.gz.tbi", chroms=CHROMS),
            expand("results/vcf/{sample}/chrStrCheck/alignmentFileHeader.out", sample=samples["sample"]),

    ##### Modules #####


    include: "rules/fq.smk"

    include: "rules/fastQC.smk"

    include: "rules/concatRawReads.smk"

    include: "rules/ref.smk"

    include: "rules/mapping.smk"

    #include: "rules/nygcAnnotations-bam.smk"

    #include: "rules/variantCalling.smk"
    include: "rules/variantCalling-Discovery.smk"

    include: "rules/dbsnpAnnotation.smk"

    include: "rules/dbsnpsSeq.smk"


**Stored sequence reads**

.. image:: img/coming-soon-storedSeq.png
   :align: center
   :width: 800

CRAM stored sequence reads - file found in ``./Iliad/workflow/cram_Snakefile``

You can delete the contents of ``./Iliad/workflow/cram_Snakefile`` and copy/paste below to replace it:
.. code-block:: console

    singularity: "library://ncherric/iliad/igdp-container:v1.16"

    include: "rules/common.smk"
    #splits=config['NYGC']['numberOfSplitRegionsFiles']

    ##### Target rules #####
    rule all:
        input: 
            expand("results/vcf/Merged-chr{chroms}-allSamples-with-rsIDs.vcf.gz.tbi", chroms=CHROMS),
            expand("results/vcf/{sample}/chrStrCheck/alignmentFileHeader.out", sample=cramSamples["cramSample"]),

    ##### Modules #####

    include: "rules/cram.smk"

    include: "rules/ref.smk"

    #include: "rules/nygcAnnotations-cram.smk"

    #include: "rules/cram_variantCalling.smk"
    include: "rules/cram_variantCalling-Discovery.smk"

    include: "rules/dbsnpsSeq.smk"

    include: "rules/dbsnpAnnotation.smk"


Model organism genome processing
================================

