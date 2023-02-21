from os import path
import re
import tempfile
from snakemake.shell import shell


# Extract arguments.
# extra = snakemake.params.get("extra", "")
extra_mpileup = snakemake.params.get("extra_mpileup", "")
extra_call = snakemake.params.get("extra_call", "")
extra_NormTrueFalse = snakemake.params.get("extra_NormTrueFalse", "")
extra_normalize = snakemake.params.get("extra_normalize", "")

vcfDir = snakemake.params.get("vcfDir", "")

reference = snakemake.input.fasta
# comment block below if you have to add specific path to referene in cram_variantCalling.smk file
# if isinstance(reference, str):
# 	reference = path.splitext(snakemake.input.ref)[0]
# else:
# 	reference = path.splitext(snakemake.input.ref[0])[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

Norm = extra_NormTrueFalse

if Norm == True:
	with tempfile.TemporaryDirectory() as tmp:
		shell(
			"bcftools +gtc2vcf"
			" -O u"
			" --bpm {snakemake.input.bpm}"
			" --egt {snakemake.input.egt}"
			" --gtcs {snakemake.input.gtc}"
			" --fasta-ref {reference} |"
			" bcftools sort -O u -T {tmp} |"
			" bcftools norm -f {reference}"
			" {extra_normalize}"
			" -O z"
			" -o {snakemake.output.vcf}"
		)

elif Norm == False:
	with tempfile.TemporaryDirectory() as tmp:
		shell(
			"bcftools +gtc2vcf"
			" -O u"
			" --bpm {snakemake.input.bpm}"
			" --egt {snakemake.input.egt}"
			" --gtcs {snakemake.input.gtc}"
			" --fasta-ref {reference} |"
			" bcftools sort -T {tmp} -O z"
			" -o {snakemake.output.vcf}"
		)