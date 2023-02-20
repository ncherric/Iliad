from os import path
from snakemake.shell import shell


# Extract arguments.
# extra = snakemake.params.get("extra", "")
wd = snakemake.params.get("workdirPath", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# shell(
#     "bcftools view"
#     " -R {snakemake.input.dbSNPrsids}"
#     " -O v"
#     " -o {snakemake.output.vcf4}"
# )

shell(
    "bcftools view"
    " --include ID==@{wd}{snakemake.input.rsidList}"
    " -O v"
    " -o {snakemake.output.vcf4}"
    " {snakemake.input.vcf3}"
)
