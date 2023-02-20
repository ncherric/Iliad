from os import path
from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get("extra", "")
extra_mpileup = snakemake.params.get("extra_mpileup", "")
extra_call = snakemake.params.get("extra_call", "")
extra_NormTrueFalse = snakemake.params.get("extra_NormTrueFalse", "")
extra_normalize = snakemake.params.get("extra_normalize", "")

sample = snakemake.params.get("sample", "")
chromsWC = snakemake.params.get("chroms_from_wc", "")
splitWC = snakemake.params.get("split_from_wc", "")

reference = snakemake.input.ref
# comment block below if you have to add specific path to referene in cram_variantCalling.smk file
# if isinstance(reference, str):
# 	reference = path.splitext(snakemake.input.ref)[0]
# else:
# 	reference = path.splitext(snakemake.input.ref[0])[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"mkdir -p results/vcf/{sample}/TempFiles/chr{chromsWC}/{splitWC}"
)

def findSequenceName(Filename,Value):
	with open(Filename) as f:
		if Value in f.read(): # read if value is in file 
			print('chr1 is present, thus, chr string exists')
			print('Adding chr to NYGC regions file')
			return True
		else:
			return False
		f.close()

Norm = extra_NormTrueFalse

filename = snakemake.input.chrStrCheck
value = 'SN:chr1'

chrStrYes = findSequenceName(filename,value)

if Norm == True:
	if chrStrYes == True:
		# keeps the chr string for "chr1"
		shell(
			"bcftools mpileup"
			" -f {reference}"
			" {extra_mpileup}"
			" -R {snakemake.input.regionsFileChr}"
			" -O u"
			" results/cram/{sample}.cram |"
			" bcftools call"
			" {extra_call}"
			" -O u |"
			" bcftools sort"
			" -m 2G"
			" -T results/vcf/{sample}/TempFiles/chr{chromsWC}/{splitWC}"
			" -O u |"
			" bcftools norm -f {reference}"
			" {extra_normalize}"
			" -O z"
			" -o {snakemake.output.splitChrVCF}"
		)

	elif chrStrYes == False:
		# removes the chr string for "1"
		shell(
			"bcftools mpileup"
			" -f {reference}"
			" {extra_mpileup}"
			" -R {snakemake.input.regionsFile}"
			" -O u"
			" results/cram/{sample}.cram |"
			" bcftools call"
			" {extra_call}"
			" -O u |"
			" bcftools sort"
			" -m 2G"
			" -T results/vcf/{sample}/TempFiles/chr{chromsWC}/{splitWC}"
			" -O u |"
			" bcftools norm -f {reference}"
			" {extra_normalize}"
			" -O z"
			" -o {snakemake.output.splitChrVCF}"
		)

elif Norm == False:
	if chrStrYes == True:
		# keeps the chr string for "chr1"
		shell(
			"bcftools mpileup"
			" -f {reference}"
			" {extra_mpileup}"
			" -R {snakemake.input.regionsFileChr}"
			" -O u"
			" results/cram/{sample}.cram |"
			" bcftools call"
			" {extra_call}"
			" -O u |"
			" bcftools sort"
			" -m 2G"
			" -T results/vcf/{sample}/TempFiles/chr{chromsWC}/{splitWC}"
			" -O z"
			" -o {snakemake.output.splitChrVCF}"
		)

	elif chrStrYes == False:
		# removes the chr string for "1"
		shell(
			"bcftools mpileup"
			" -f {reference}"
			" {extra_mpileup}"
			" -R {snakemake.input.regionsFile}"
			" -O u"
			" results/cram/{sample}.cram |"
			" bcftools call"
			" {extra_call}"
			" -O u |"
			" bcftools sort"
			" -m 2G"
			" -T results/vcf/{sample}/TempFiles/chr{chromsWC}/{splitWC}"
			" -O z"
			" -o {snakemake.output.splitChrVCF}"
		)

