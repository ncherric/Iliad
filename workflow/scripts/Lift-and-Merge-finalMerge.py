from os import path
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "bcftools merge -l {snakemake.input.mergeList}"
    " -m both"
    " --force-samples"
    " -O z"
    " -o {snakemake.output.mergeVCF}"
)

