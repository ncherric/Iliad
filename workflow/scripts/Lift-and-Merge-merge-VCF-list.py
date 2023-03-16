from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")
project = snakemake.params.get("project_from_wc")
refAssemblyVersion = snakemake.params.get("refAssemblyVersion_from_wc")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

sed_arg = [ "sed 's/ /\\n/g'" ]

shell(
	"echo {snakemake.input.VCFs_in_37} |"
	" {sed_arg} - > data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/VCFlist-from-37path.txt"
	)

shell(
	"echo {snakemake.input.VCFs_in_38} |"
	" {sed_arg} - > data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/VCFlist-from-38path.txt"
	)

shell(
	"cat data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/VCFlist-from-37path.txt data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/VCFlist-from-38path.txt "
	"> {snakemake.output.mergeList}"
)

