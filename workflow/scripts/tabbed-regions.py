from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
chrStr=str("chr")
blank=str('""')

# removes the chr string for "1"
shell(
    "sed 's/ /\t/g' {snakemake.input.uniqFile} |"
    " sed '1d' - |"
    " awk '{{text=$1; gsub(/chr/, {blank}); print }}' OFS='\t' -"
    " > {snakemake.output.tabFile}"
)

