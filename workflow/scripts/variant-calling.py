from os import path
from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get("extra", "")

reference = snakemake.input.ref
if isinstance(reference, str):
	reference = path.splitext(snakemake.input.ref)[0]
else:
	reference = path.splitext(snakemake.input.ref[0])[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)



shell(
    "bcftools mpileup -d 8000"
    " -f {reference}"
    " -B"
    " -R {snakemake.input.regionsFile}"
    " -O u"
    " {snakemake.input.sortedbam} |"
    " bcftools call"
    " -m"
    " -A"
    " -O z"
    " -o {snakemake.output.splitChrVCF}"
)
