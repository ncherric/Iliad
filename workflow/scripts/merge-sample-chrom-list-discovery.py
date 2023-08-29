from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
sed_arg = [ "sed 's/ /\\n/g'" ]

shell(
    "echo {snakemake.input.ChrVCF} |"
    " {sed_arg} - > {snakemake.output.mergeList}"
)

