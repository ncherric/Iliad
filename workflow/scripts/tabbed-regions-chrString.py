from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# keeps the chr string for "chr1"
shell(
    "sed 's/ /\t/g' {snakemake.input.uniqFile} |"
    " sed '1d' -"
    " > {snakemake.output.tabChrFile}"
)

