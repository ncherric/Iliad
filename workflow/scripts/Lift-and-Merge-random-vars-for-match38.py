from os import path
from snakemake.shell import shell
import pandas as pd

# Extract arguments.
project=snakemake.params.get("project_from_wc")
refAssemblyVersion=snakemake.params.get("refAssemblyVersion_from_wc")
vcf38=snakemake.params.get("vcf38_from_wc")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"mkdir -p data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/"
)

shell(
	"bcftools view -O v -o {snakemake.output.vcfForRandomizing} {snakemake.input.file38[0]}"
)
shell(
	"awk '/^[^#]/ {{print $1,$2,$3}}' OFS='\t' {snakemake.output.vcfForRandomizing} |"
	" grep -E 'rs' - |"
	" shuf -n 3 - >"
	" {snakemake.output.myData_randomVarsForMatch}"
)
