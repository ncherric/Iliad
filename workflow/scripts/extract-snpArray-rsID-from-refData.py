from os import path
from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "bcftools view"
    " -R {snakemake.input.keeperSNPs}"
    " -O z"
    " -o {snakemake.output.RefChrVCFoverlap}"
    " {snakemake.input.RefChrVCF}"
)