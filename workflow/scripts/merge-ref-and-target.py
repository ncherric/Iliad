from os import path
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"bcftools merge {snakemake.input.ReferenceVCF} {snakemake.input.TargetVCF}"
	" -m both"
	" -O z"
	" -o {snakemake.output.MergedTargetRef}"
)
