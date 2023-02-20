#wget text file from illumina website
#unzip
#convert text to JSON
#Import pyVCF and then import vcf.gz to replace the loci names in the ID column

import wget
import zipfile
import re
from os import path, chmod
import json
from datetime import datetime
import csv
import pysam

import re
import tempfile
from snakemake.shell import shell


# Extract arguments.
# extra = snakemake.params.get("extra", "")

url = snakemake.params.get("rsidsUrl", "")
rzip = snakemake.params.get("rzip", "")
workdir = snakemake.params.get("workdirPath", "")
rfile = snakemake.params.get("rfile", "")

location = workdir + 'illumina/support_files/' + rzip
rloc = workdir + 'illumina/support_files/'
rfileLoc = workdir + 'illumina/support_files/' + rfile
vcfLoc = workdir + 'data/snp_array/vcf/'

# 1 - IF import wget and import pysam would work on snakemake.....
# Download the conversion file
print('Beginning Loci Name to rsID' + '/n' + 'Conversion file download from illumina with wget module')
#url = 'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/multiethnic-global/multi-ethnic-global-8-d1-rsids.zip'
# wget.download(url, abspath(join(dirname(__file__), rzip)))
wget.download(url, out = location)

# unzip the conversion file
with zipfile.ZipFile(location, 'r') as zip_ref:
	conversionFile = zip_ref.namelist()
	print(f'\n' + f'filename is: {conversionFile}')
	chmod(rloc + str(conversionFile)[2:-2], 0o777)
	zip_ref.extractall(rloc)
conversionFile = str(conversionFile)[2:-2]
print(f'\n' + f'filename as string: {conversionFile}')

# # 2
# # convert loci name to rsid conversion text file to JSON file
dict_conversion = {}
# creating dictionary
with open(rfileLoc) as cf:
	for i in cf:
		# reads each line and trims off extra spaces  
		# and gives only the valid words 
		command, description = i.strip().split(None, 1)
		dict_conversion[command] = description.strip()

# creating json file
# the JSON file is named as test1 
json_file = open(vcfLoc + "A-lociname-to-rsid-conversion.json", "w")
json.dump(dict_conversion, json_file, indent = 4, sort_keys = False) 
json_file.close()

# removing key == value pairs in dictionary
key_num1 = len(dict_conversion)
print(f'Keys count BEFORE removal of key==value instances: {key_num1}')
for key in list(dict_conversion):
	if key == dict_conversion[key]:
		del dict_conversion[key]
key_num2 = len(dict_conversion)
print(f'Keys count AFTER removal of key==value instances: {key_num2}')


# opened the vcf file and the new file to write via csv library to one column text file
###
### Needs more dev on how to output to vcf.gz file via pysam or output to vcf via reading in line by line of vcf

vcf_in = pysam.VariantFile(snakemake.input.vcf2)  # auto-detect input format
#vcf_out = VariantFile("test/out-test.vcf.gz", 'w', header=vcf_in.header)

lociNames = []
with open(snakemake.output.rsidsConversionFile, 'w') as new_file:
	count = 0
	for rec in vcf_in.fetch():
		count += 1
		#print(f'rec id line = {rec.id}')
		#print('type is: {}'.format(type(rec.id)))
		lociNames.append(rec.id)
		if(count % 10000 == 0):
			print('lociNames append count is: {} {}'.format(count, datetime.now()))
	count = 0
	for loci in lociNames:
		count += 1
		if loci in dict_conversion:
			#print('type = {}, lociName = {}, replace rs-id = {}'.format(type(loci), loci, dict_conversion[loci]))
			new_file.write(dict_conversion[loci]+'\n')
		else:
			#print('line not found')
			new_file.write(loci+'\n')
		if(count % 100000 == 0):
			print('{} {}'.format(count, datetime.now()))

print('done')
