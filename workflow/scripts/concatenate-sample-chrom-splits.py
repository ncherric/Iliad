from os import path
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")
# TempDir = snakemake.params.get("sampleChrTempDir", "")
sampleWC = snakemake.params.get("samples_from_wc", "")
chromsWC = snakemake.params.get("chroms_from_wc", "")
# indexConcatOutput = snakemake.params.get("indexConcatOutput", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"mkdir -p results/vcf/{sampleWC}/TempFiles/chr{chromsWC}/"
)

shell(
	"bcftools concat -O u"
	" -f {snakemake.input.concatList} |"
	" bcftools sort"
	" -m 8G"
	" -T results/vcf/{sampleWC}/TempFiles/chr{chromsWC}/"
	" -O z"
	" -o {snakemake.output.concatChrVCF}"
)

