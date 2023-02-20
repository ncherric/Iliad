from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "bcftools view  --include ID!=@{snakemake.input.GenTrain_ToDropSNPs}"
    " -O u {snakemake.input.vcf4} |"
    " bcftools view --include ID!=@{snakemake.input.ClusterSep_ToDropSNPs}"
    " -O v -o {snakemake.output.vcf5}"
)

shell(
    "bcftools view  --include ID!=@{snakemake.input.GenTrain_ToDropSNPs}"
    " -O u {snakemake.input.vcf4} |"
    " bcftools view --include ID!=@{snakemake.input.ClusterSep_ToDropSNPs}"
    " -O z -o {snakemake.output.vcf5gz}"
)
