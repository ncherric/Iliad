from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "bcftools query"
    " -f '%CHROM\t%POS\t%REF\t\n'"
    " {snakemake.input.startVCF} |"
    " shuf -n 5"
    " > {snakemake.output}"
)
