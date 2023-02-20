from os import path
from snakemake.shell import shell
import pandas as pd
from builtins import str
from io import open


# Extract arguments.
project = snakemake.params.get("project_from_wc")
refAssemblyVersion = snakemake.params.get("refAssemblyVersion_from_wc")

versionCheck = snakemake.input.versionCheck
dbSNP_Extracted = snakemake.input.dbSNP_Extracted
dbSNP_Matches = snakemake.output.dbSNP_Matches
VCF38 = snakemake.params.get("VCF38", "")


shell(
	"touch {snakemake.output.tmpoFile38path}"
)

# versionCheckEnc = versionCheck.encode()
# dbSNP_ExtractedEnc = dbSNP_Extracted.encode()
# dbSNP_MatchesEnc = dbSNP_Matches.encode()
# , encoding="latin-1"

def findVersion(Fname, dbsnpfile, output_file):
    with open(dbsnpfile,'r') as read1:
      f1 = read1.read().split('\n')

    with open(Fname,'r') as read2:
      f2 = read2.read().split('\n')

    with open(output_file,'w', encoding="latin-1") as write1:
      for line in f2:
        if line in f1:
          write1.write(line + '\n')
        else:
          write1.write('blank\n')

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

findVersion(versionCheck,dbSNP_Extracted,dbSNP_Matches)

dfVCF = pd.read_csv(versionCheck, sep = '\t', header=None)
dfDBSNP = pd.read_csv(dbSNP_Matches, sep = '\t', header=None)
value = 'blank'

blank_InChecks = AreThere_blank_InChecks(dbSNP_Matches,value)

if blank_InChecks is True:
    print("Atleast 1 out of 3 lines in dbSNP matches file returned a blank")
    print("This was a double check validation for if file was version 38.\n")
    print("You will have to check the version of this file.\n")

    shell(
        "bcftools view -O z -o {VCF38} {snakemake.input.keep23}"
    )
    shell(
        "bcftools index --tbi {VCF38}"
    )

elif blank_InChecks is False:
    print("\nThere were not any blanks in dbSNP matches file\n")
    print("Assessed VCF input will be kept as 38, processed, and merged with other files.\n")
    df1 = dfVCF.sort_values(dfVCF.columns[2])
    df2 = dfDBSNP.sort_values(dfDBSNP.columns[2])
    df37=df2.equals(df1)

    if df37 is True:
        print("Double checked the dbSNP matches file and VCF randoms\n")
        print("They are equal and passing on 37 VCF to next step\n")
        shell(
            "bcftools view -O z -o {VCF38} {snakemake.input.keep23}"
        )
        shell(
            "bcftools index --tbi {VCF38}"
        )

    elif df37 is False:
        print("Double checked the dbSNP matches file and VCF randoms\n")
        print("There was some discrepancy in the two files, but not a blank\n")
        exit("Exiting the program. You will have to check your files!")
