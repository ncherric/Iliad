from os import path
from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get("extra", "")
sample = snakemake.params.get("sample", "")

chromsWC = snakemake.params.get("chroms_from_wc", "")
splitWC = snakemake.params.get("split_from_wc", "")

reference = snakemake.input.ref
# if isinstance(reference, str):
# 	reference = path.splitext(snakemake.input.ref)[0]
# else:
# 	reference = path.splitext(snakemake.input.ref[0])[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "mkdir -p results/vcf/{sample}/TempFiles/chr{chromsWC}/{splitWC}"
)


shell(
    "bcftools mpileup -d 8000"
    " -f {reference}"
    " -B"
    " -R {snakemake.input.regionsFile}"
    " -O u"
    " results/downloads/{sample}/{sample}.final.cram |" # remove .final and add a renaming step to download to soft code download name
    " bcftools call"
    " -m"
    " -A"
    " -O u |"
    " bcftools sort"
    " -m 2G"
    " -T results/vcf/{sample}/TempFiles/chr{chromsWC}/{splitWC}"
    " -O z"
    " -o {snakemake.output.splitChrVCF}"
)
