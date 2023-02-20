from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

reference = snakemake.input.fasta
if isinstance(reference, str):
	reference = path.splitext(snakemake.input.fasta)[0]
else:
	reference = path.splitext(snakemake.input.fasta[0])[0]


log = snakemake.log_fmt_shell(stdout=False, stderr=True)


with tempfile.TemporaryDirectory() as tmp:
	shell(
		"bcftools view "
		"-m2 -M2 -v snps "
		"-R {snakemake.input.regionsFile} "
		"-O u "
		"{snakemake.input.vcf} | "
		"bcftools sort "
		"-m 3G "
		"-T {tmp} "
		"-O u | "
		"bcftools norm "
		"-c w "
		"-f {reference} "
		"-O u | "
		# "bcftools norm "
		# "-m +both "
		# "-O u | "
		"bcftools norm -d both "
		"-O v "
		"-o {snakemake.output.vcf2}"
	)
