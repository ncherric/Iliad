from os import path
from snakemake.shell import shell
import pandas as pd

# Extract arguments.
project=snakemake.params.get("project_from_wc")
refAssemblyVersion=snakemake.params.get("refAssemblyVersion_from_wc")
vcf=snakemake.params.get("vcf_from_wc")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"mkdir -p data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/"
)

shell(
	"bcftools view -O v -o {snakemake.output.vcfForRandomizing} {snakemake.input.keep23}"
)
shell(
	"awk '/^[^#]/ {{print $1,$2,$3}}' OFS='\t' {snakemake.output.vcfForRandomizing} |"
	" grep -E 'rs' - |"
	" shuf -n 3 - >"
	" {snakemake.output.myData_randomVarsForMatch}"
)
