# TL;DR setup


## General Input and Output

| Input                                 | Output                                     |
|:-------------------------------------:|:------------------------------------------:|
| FASTQ data or FTP links to FASTQ data | quality-controlled VCF for each chromosome |


## **Please make sure that your conda environment for Iliad is activated** - ``conda activate iliadEnv`` or ``mamba activate iliadEnv``

**Modify the configuration file** ``workdirPath`` parameter to the appropriate path leading up to and including ``/Iliad`` and a final forward slash e.g. ``/Path/To/Iliad/``. 
The configuration file is found in ``config/config.yaml``.

```yaml

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
```

### You might consider changing some other parameters to your project needs that are pre-set and include:

Homo sapiens GRCh38 release 104 reference genome

```yaml

    ref:
      species: homo_sapiens
      release: 104
      build: GRCh38
```

Use an Excel sheet or CSV file with no header and the following two columns/fields:

```console

    Sample   Unique sample identifier
    URL   raw sequence data download FTP link
```

#### Example: **UserSampleTable.xlsx** or **UserSampleTable.csv** are found in the ``/Iliad/config/`` directory

| KPGP-00127 | ftp://ftp.kobic.re.kr/pub/KPGP/2020_release_candidate/WGS_SR/KPGP-00127/KPGP-00127_L1_R1.fq.gz |
| --- | --- |
| KPGP-00127 | ftp://ftp.kobic.re.kr/pub/KPGP/2020_release_candidate/WGS_SR/KPGP-00127/KPGP-00127_L1_R2.fq.gz |

This exact template exists already in ``/Iliad/config/UserSampleTable.xlsx``. (The Excel Viewer extension on VS code is really handy for editing!)
If you already have the sequence files and are not downloading open-source data, you have the option to place your data into the ``Iliad/data/fastq/`` directory.

Whether you are automatically downloading via Iliad or you manually place data into ``Iliad/data/fastq/`` directory,
you need to provide a separate ``samples.tsv`` file where the TSV file has a header line with only one field named ``sample``.

```console

    sample  HEADER
    SAMPLE1 sample identifier
    SAMPLE2 sample identifier
```

#### Example: **samples.tsv** found in the ``/Iliad/config/`` directory

| sample |
| --- |
| KPGP-00127 |

since this module is the main snakefile, Snakemake will automatically detect it without the flag. 
(Please make sure that your conda environment for Iliad is activated - ``conda activate iliadEnv`` or ``mamba activate iliadEnv``)

```console

    $ snakemake --cores 1
```

and combined with other user-specified snakemake flags such as ``--cores``.

If you plan to use on a local machine or self-built server without a job scheduler the default command to run is the following:

```console

   $ snakemake -p --use-singularity --use-conda --cores 1 --jobs 1 --default-resource=mem_mb=10000 --latency-wait 120
```

However, there is a file included in the ``Iliad`` directory named - ``snakemake.sh`` that will be useful in batch job submission. 
Below is an example snakemake workflow submission in SLURM job scheduler. 
Please read the shell variables at the top of the script and customize to your own paths and resource needs.

```console

   $ sbatch snakemake.sh
```

If you would like more in-depth information and descriptions, please continue to the next sections below. 
Otherwise, you have completed the TL;DR setup section.