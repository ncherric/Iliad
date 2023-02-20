from os import path
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"mkdir -p data/target_ref_merge/TempFiles/
)

shell(
	"bcftools concat -O u"
	" -f {snakemake.input.concatList} |"
	" bcftools sort"
	" -m 8G"
	" -T data/target_ref_merge/TempFiles/"
	" -O z"
	" -o {snakemake.output.concatVCF}"
)

