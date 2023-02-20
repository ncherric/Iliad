from builtins import str
import os
from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


# could add smaller functions for sub processes happening in conditional statements

# def autoDetect_Zip_Unzip(vcfmergedir, liftoverbool, version38, version37, dbsnpDir, dbsnpFile):
	# Checking each file type extension and performing necessary tasks per .vcf or .vcf.gz
	# iterate over files in that directory
# folder = MergeDir
chrStr=str('"chr"')
blank=str('""')
redefineHeader=str('"contig=<ID="')

print("The file - {} - is a VCF\n".format(snakemake.input.vcfInput))

shell(
    "awk '{{ gsub(/chr/,{blank}, $1); print }}' OFS='\t' {snakemake.input.vcfInput} |"
    " awk '{{text=$1; gsub(/contig=<ID=chr/,{redefineHeader}); print }}' OFS='\t'"
    " > {snakemake.output.NoChrString}"
)