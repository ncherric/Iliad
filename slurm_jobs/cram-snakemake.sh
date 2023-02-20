#! /bin/bash
# ADD YOUR personal mail-user, output log file, and job-name
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=1-12:30:0
#SBATCH --mail-user=
#SBATCH --mail-type=ALL,TIME_LIMIT
#SBATCH --output=
#SBATCH --job-name=

# Make sure singularity is accessible in your HPC environment - i.e. module load singularity

# Change into your Working Iliad Directory - uncomment below
# cd 

# If you receive an error about conda init, you may have to initialize conda for your bash or sh shell - uncomment below
# conda init bash

# you want one process that can use 16 cores for multithreading: --ntasks=1 --cpus-per-task=16
sbcmd="sbatch --ntasks=1 --cpus-per-task={threads} --mem={resources.mem_mb}"
sbcmd+=" --time=1-12:0:0 --output=[Working Iliad Directory]/logs/{rule}.{wildcards}.o"
sbcmd+=" --error=[Working Iliad Directory]/Iliad/logs/{rule}.{wildcards}.e"
sbcmd+=" --mail-user= --mail-type=ALL,TIME_LIMIT"

snakemake -p --use-singularity --use-conda --cores 8 --jobs 8 --snakefile workflow/cram_Snakefile --cluster "$sbcmd" --latency-wait 60