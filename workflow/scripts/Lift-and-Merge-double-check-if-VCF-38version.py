from os import path
from snakemake.shell import shell
import pandas as pd

# Extract arguments.
project=snakemake.params.get("project_from_wc")
refAssemblyVersion=snakemake.params.get("refAssemblyVersion_from_wc")
dbsnpDir=snakemake.params.get("dbsnpDir", "")
dbsnpFile=snakemake.params.get("dbsnpFile", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"mkdir -p data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/"
)

shell(
	"zcat {snakemake.input.keep23} |"
	" awk '/^[^#]/ {{print $1,$2,$3}}' OFS='\t' - |"
	" shuf -n 3 - >"
	" {snakemake.output.versionCheck}"
)
# now i have a file with 3 random "chr	pos"

# need to get those 3 from dbSNP file
shell(
	"bcftools view -R {snakemake.output.versionCheck}"
	" -O v"
	" -o data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/tmp.file.txt"
	" {dbsnpDir}{dbsnpFile}"
)

shell(
	" awk '/^[^#]/ {{print $1,$2,$3}}' OFS='\t' data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/tmp.file.txt >"
	" {snakemake.output.dbSNP_Extracted}"
)
