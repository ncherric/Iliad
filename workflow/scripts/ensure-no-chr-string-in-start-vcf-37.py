from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "awk '{ gsub(/chr/,"", $1); print }' OFS='\t' {snakemake.input[0]} "
    " > {snakemake.output.startVCFwChr}"
)











        
