from os import path
from snakemake.shell import shell
import pandas as pd
from builtins import str
from io import open


# Extract arguments.
project = snakemake.params.get("project_from_wc")
refAssemblyVersion = snakemake.params.get("refAssemblyVersion_from_wc")
vcf = snakemake.params.get("vcf_from_wc")

myData_randomVarsForMatch = snakemake.input.myData_randomVarsForMatch
dbSNP_ExtractedVarsForMatch = snakemake.input.dbSNP_ExtractedVarsForMatch
dbSNP_Matches = snakemake.input.dbSNP_Matches

def AreThere_blank_InChecks(Filename, val):
	with open(Filename,'r') as f:
		if val in f.read():
			print('\nAtleast 1/3 SNPs was blank')
			print('Cross referencing with 38 dbsnp now')
			return True
		else:
			print('\nAll 3 SNPs were there, this is most likely GRCh37')
			print('Proceeding to annotate, lift over to 38, and merge')
			return False


dfVCF = pd.read_csv(myData_randomVarsForMatch, sep = '\t', header=None)
dfDBSNP = pd.read_csv(dbSNP_Matches, sep = '\t', header=None)
value = 'blank'

blank_InChecks = AreThere_blank_InChecks(dbSNP_Matches,value)

# this is the 38 destination path
if blank_InChecks is True:
	print("Atleast 1 out of 3 lines in dbSNP matches file returned a blank - a non-match or missing extracted site")
	print("Assessed VCF input will be further assessed with 38 dbSNP checkpoint\n")
# making the 37 filename and basename files (these will be empty for v38s)
	shell(
		"touch {snakemake.output.fileWithPathIfVersion37}"
	)

	shell(
		"echo > {snakemake.output.fileWithBaseNameIfVersion37}"
	)

# making the 38 filename and basename files

	shell(
		"echo '38path_{vcf}' > {snakemake.output.fileWithBaseNameIfVersion38}"
	)

	shell(
		"echo 'data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/38path_{vcf}.vcf.gz' > {snakemake.output.fileWithPathIfVersion38}"
	)


# make 38 directory

	shell(
		"mkdir -p data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/"
	)


# making the important next step VCFs for 38
	shell(
		"bcftools view -O z -o data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/38path_{vcf}.vcf.gz {snakemake.input.keep23}"
	)

	shell(
		"bcftools index --tbi data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/38path_{vcf}.vcf.gz"
	)

# this is the 37 destination path
elif blank_InChecks is False:
	print("\nThere were not any blanks in dbSNP matches file\n")
	print("Assessed VCF input will be converted from 37 to 38\n")
	df1 = dfVCF.sort_values(dfVCF.columns[2])
	df2 = dfDBSNP.sort_values(dfDBSNP.columns[2])
	df37=df2.equals(df1)

	if df37 is True:
		print("Double checked the dbSNP matches file and VCF randoms\n")
		print("They are equal and passing on 37 VCF to next step\n")

# making the 38 filename and basename files(these will be empty for v37s)
		shell(
			"touch {snakemake.output.fileWithPathIfVersion38}"
		)

		shell(
			"echo > {snakemake.output.fileWithBaseNameIfVersion38}"
		)

# making the 37 filename and basename files

		shell(
			"echo '37path_{vcf}' > {snakemake.output.fileWithBaseNameIfVersion37}"
		)

		shell(
			"echo 'data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3A-InputVCFs-37/37path_{vcf}.vcf.gz' > {snakemake.output.fileWithPathIfVersion37}"
		)

# make 37 directory

		shell(
			"mkdir -p data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3A-InputVCFs-37/"
		)


# making the important next step VCFs for 37
		shell(
			"bcftools view -O z -o data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3A-InputVCFs-37/37path_{vcf}.vcf.gz {snakemake.input.keep23}"
		)

		shell(
			"bcftools index --tbi data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3A-InputVCFs-37/37path_{vcf}.vcf.gz"
		)

	elif df37 is False:
		print("Double checked the dbSNP matches file and VCF randoms\n")
		print("There was some discrepancy in the two files, but not a blank\n")
		exit("Exiting the program. You will have to check your files!")

