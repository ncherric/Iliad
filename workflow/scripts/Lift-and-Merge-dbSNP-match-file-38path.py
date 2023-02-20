from os import path
from snakemake.shell import shell
import pandas as pd
from builtins import str
from io import open


# Extract arguments.

myData_randomVarsForMatch = snakemake.input.myData_randomVarsForMatch
dbSNP_ExtractedVarsForMatch = snakemake.input.dbSNP_ExtractedVarsForMatch
dbSNP_Matches = snakemake.output.dbSNP_Matches

def findVersion(Fname, dbsnpfile, output_file):
    with open(dbsnpfile,'r', encoding="latin-1") as read1:
      f1 = read1.read().split('\n')

    with open(Fname,'r', encoding="latin-1") as read2:
      f2 = read2.read().split('\n')

    with open(output_file,'w', encoding="latin-1") as write1:
      for line in f2:
        if line in f1:
          write1.write(line + '\n')
        else:
          write1.write('blank\n')

findVersion(myData_randomVarsForMatch,dbSNP_ExtractedVarsForMatch,dbSNP_Matches)
