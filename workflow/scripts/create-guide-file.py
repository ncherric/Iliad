from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
chrStr=str("chr")
blank=str('""')

shell(
    "awk '{{print $1,$2,$3}}' {snakemake.input.dbsnpExtractedIDsFile} |"
    " awk '{{text=$1; gsub(/chr/, {blank}); print }}' OFS='\t' - |"
    " awk '$1={chrStr}$1' OFS='\t'"
    " > {snakemake.output.guideFile}"
)

    # " awk '{{ gsub(/chr/, "", $1); print }}' OFS='\t' - |"
