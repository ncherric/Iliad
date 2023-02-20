from os import path
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")
# indexOutput = snakemake.params.get("indexOutput", "")

dbsnpFile = snakemake.params.get("dbsnpFile", "")
dbsnpDir = snakemake.params.get("dbsnpDir", "")
workdirPath = snakemake.params.get("workdirPath", "")
# tmpoFile = snakemake.input.tmpoFile

dbSNP = str(workdirPath + dbsnpDir + dbsnpFile)

# rmFile = str(workdirPath + {snakemake.input.tmpoFile})
# rmFile = str(workdirPath + tmpoFile)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "bcftools annotate"
    " -a {dbSNP}"
    " -c ID"
    " -O z"
    " -o {snakemake.output.mergeChrVCF_rsID}"
    " {snakemake.input.mergeChrVCF}"
    # " && bcftools index --tbi {snakemake.output.mergeChrVCF_rsID}"
    # " && rm -f {rmFile}"
)
