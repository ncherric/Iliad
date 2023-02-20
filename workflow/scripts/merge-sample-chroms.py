from os import path
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")
indexMergeOutput = snakemake.params.get("indexMergeOutput", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "bcftools merge -l {snakemake.input.mergeList}"
    " -m both"
    " -O z"
    " -o {snakemake.output.mergeChrVCF}"
    # " && sleep 15s"
    # " && bcftools index --tbi {snakemake.output.mergeChrVCF}"
)

