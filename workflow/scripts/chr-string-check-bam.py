from os import path
from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get("extra", "")
# indexOutput = snakemake.params.get("indexOutput", "")
sample = snakemake.params.get("sample", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "mkdir -p results/vcf/{sample}/chrStrCheck/"
)

shell(
	"/usr/bin/samtools-1.15/bin/samtools head {snakemake.input.sortedBam} |"
	" awk '{{print $1,$2}}' - |"
	" grep -w '@SQ' - >"
	" results/vcf/{sample}/chrStrCheck/alignmentFileHeader.out"
)


