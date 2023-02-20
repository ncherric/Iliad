from os import path
from snakemake.shell import shell
import pandas as pd
from builtins import str
from io import open


# Extract arguments.
project = snakemake.params.get("project_from_wc")
refAssemblyVersion = snakemake.params.get("refAssemblyVersion_from_wc")
vcf38 = snakemake.params.get("vcf38_from_wc")

myData_randomVarsForMatch = snakemake.input.myData_randomVarsForMatch
dbSNP_ExtractedVarsForMatch = snakemake.input.dbSNP_ExtractedVarsForMatch
dbSNP_Matches = snakemake.input.dbSNP_Matches

def AreThere_blank_InChecks(Filename, val):
	with open(Filename,'r') as f:
		if val in f.read():
			print('\nAtleast 1/3 SNPs was blank')
			print('the double validation failed')
			return True
		else:
			print('\nAll 3 SNPs were there, this is most likely GRCh38')
			print('Proceeding to annotate as version 38, and merge with others')
			return False


dfVCF = pd.read_csv(myData_randomVarsForMatch, sep = '\t', header=None)
dfDBSNP = pd.read_csv(dbSNP_Matches, sep = '\t', header=None)
value = 'blank'

blank_InChecks = AreThere_blank_InChecks(dbSNP_Matches,value)

# this is the errored out, non validated destination path
if blank_InChecks is True:
	print("Atleast 1 out of 3 lines in dbSNP matches file returned a blank - a non-match or missing extracted site")
	print("Assessed VCF input will be further removed from this workflow. Failed second validation\n")
# making the not validated filename and basename files (these will be empty for v38s)
	shell(
		"touch {snakemake.output.fileWithPathIfVersionNotValidated}"
	)

	shell(
		"echo > {snakemake.output.fileWithBaseNameIfVersionNotValidated}"
	)

# making the 38 filename and basename files

	shell(
		"echo > {snakemake.output.fileWithBaseNameIfVersion38}"
	)

	shell(
		"echo > {snakemake.output.fileWithPathIfVersion38}"
	)


# this is the validated 38 destination path
elif blank_InChecks is False:
	print("\nThere were not any blanks in dbSNP matches file\n")
	print("Assessed VCF input will be processed in the 38 path\n")
	df1 = dfVCF.sort_values(dfVCF.columns[2])
	df2 = dfDBSNP.sort_values(dfDBSNP.columns[2])
	df38=df2.equals(df1)

	if df38 is True:
		print("Double checked the dbSNP matches file and VCF randoms\n")
		print("They are equal and passing on 37 VCF to next step\n")

# making the errored out filename and basename files(these will be empty for v38s)
		shell(
			"touch {snakemake.output.fileWithPathIfVersionNotValidated}"
		)

		shell(
			"echo > {snakemake.output.fileWithBaseNameIfVersionNotValidated}"
		)

# making the 38 filename and basename files

		shell(
			"echo 'valid_{vcf38}' > {snakemake.output.fileWithBaseNameIfVersion38}"
		)

		shell(
			"echo 'data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/validated/valid_{vcf38}.vcf.gz' > {snakemake.output.fileWithPathIfVersion38}"
		)

# make 37 directory

		shell(
			"mkdir -p data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/validated/"
		)


# making the important next step VCFs for 37
		shell(
			"bcftools view -O z -o data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/validated/valid_{vcf38}.vcf.gz {snakemake.input.file38}"
		)

		shell(
			"bcftools index --tbi data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/validated/valid_{vcf38}.vcf.gz"
		)

	elif df38 is False:
		print("Double checked the dbSNP matches file and VCF randoms\n")
		print("There was some discrepancy in the two files, but not a blank\n")
		exit("Exiting the program. You will have to check your files!")

