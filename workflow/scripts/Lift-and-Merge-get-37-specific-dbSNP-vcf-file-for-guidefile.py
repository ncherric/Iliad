from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

dbsnpDir=snakemake.params.get("dbsnpDir", "")
dbsnpFile=snakemake.params.get("dbsnpFile", "")

shell(
	"bcftools view -i ID=@{snakemake.input.combinedSNPlist}"
	" -O v"
	" -o {snakemake.output.projectSpecificDBSNPvcf}"
	" {dbsnpDir}{dbsnpFile}"
)