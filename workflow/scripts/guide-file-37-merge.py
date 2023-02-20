from builtins import str
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

chrStr=str('"chr"')
blank=str('""')

shell(
    "awk '{{print $1,$2,$3}}' {snakemake.input.dbsnpExtractedIDsFile} |"
    " awk '{{ gsub(/chr/,{blank}, $1); print }}' OFS='\t' -"
    " > {snakemake.output.guideFile}"
)