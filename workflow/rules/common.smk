import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
import pathlib
import os
from os.path import join
from glob import glob


###### Config file and sample sheets #####
configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

cramSamples = pd.read_table(config["cramSamples"]).set_index("cramSample", drop=False)
validate(cramSamples, schema="../schemas/cramSamples.schema.yaml")

# for merging vcf and vcf.gz - submodule
vcfs = pd.read_csv(config["vcfs"]).set_index("baseFileName_VCF", drop=False)


# add condition to run all by default or specific chr by user choice in config file
CHROMS=list(range(1, 23))
CHROMS.extend(['X'])

NAMES=config['LiftoverSub']['filename']

###### Splash Screen and ensure correct MODULE/WORKFLOW is being used #####

print(
# "__author__ = Noah Herrick   \n"
# "__copyright__ = Copyright 2023, Noah Herrick   \n"
# "__email__ = ncherric at iu dot edu   \n"
# "__license__ = MIT   \n\n"
"                                                      \n"
"              ░░░░░░░░░░░░░░░░░░░░░░░░░░░░            \n"
"          ░░░░░░▒▒▓▓▒▒▒▒▒▒░░▒▒▒▒▓▓▒▒░░░░░░░░░░        \n"
"          ░░▒▒██▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░▒▒▓▓░░░░░░░░      \n"
"    ░░░░▒▒▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░▓▓░░░░░░    \n"
"    ░░▒▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░▒▒░░░░░░  \n"
"  ░░░░▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░░░\n"
"  ░░░░██▓▓▓▓▓▓▓▓▓▓▒▒▓▓▒▒░░░░▒▒▓▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▓▓░░░░░░\n"
"   ░▒▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒░░      ░░▒▒▒▒▒▒▒▒▒▒▒▒▒▒▓▓▒▒░░░ \n"
"    ▒▒▓▓▓▓▓▓▓▓▓▓▒▒██░░▒▒▒▒▒▒  ▒▒▓▓▒▒▒▒▒▒▒▒▒▒▒▒▓▓▒▒░░  \n"
"     ▒▓▓▓▓▓▓▓▓▓▓▓▓██  ▓▓▒▒▒▒░░▒▒▓▓▓▓▒▒▓▓▒▒▒▒▒▒▓▓▒▒░   \n"
"    ▒▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒    ░░▒▒▓▓▓▓▓▓▒▒▒▒▒▒▒▒▓▓▒▒░░  \n"
"   ░▒▒██▓▓▓▓▓▓▒▒▒▒▒▒▓▓░░    ▒▒▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▓▓▒▒░░░ \n"
"  ░░░░▓▓▒▒▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒▒▒▒▒▓▓▓▓▒▒▒▒▒▒▓▓▒▒░░░░░░\n"
"    ░░▒▒▓▓▒▒▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒▒▒▓▓▓▓▓▓▒▒▒▒▓▓▒▒░░░░  \n"
"    ░░▒▒██▒▒▓▓▓▓▓▓▓▓▒▒▓▓▒▒▒▒▒▒▓▓▓▓▓▓▒▒▒▒▒▒▒▒▒▒▒▒░░░░  \n"
"      ░░░░▒▒██▒▒▓▓▓▓▓▓▓▓▓▓▓▓▒▒▓▓▓▓▒▒▒▒▒▒██▒▒▒▒░░░░    \n"
"        ░░░░▒▒▓▓▓▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▓▓██▒▒▒▒▒▒░░░░      \n"
"        ░░░░▒▒▒▒▓▓██▓▓▒▒▒▒▒▒▒▒▓▓██▓▓▒▒▒▒▒▒░░░░        \n"
"          ░░░░▒▒▒▒▒▒▓▓▓▓▓▓██▓▓▓▓▒▒▒▒▒▒░░░░░░          \n"
"              ▒▒░░░░░░░░░░▒▒▒▒▒▒░░░░░░░░              \n\n"

"         ██╗  ██╗       ██╗   █████╗   ██████╗ \n"
"         ██║  ██║       ██║  ██╔══██╗  ██╔══██╗\n"
"         ██║  ██║       ██║  ███████║  ██║  ██║\n"
"         ██║  ██║       ██║  ██╔══██║  ██║  ██║\n"
"         ██║  ███████╗  ██║  ██║  ██║  ██████╔╝\n"
"         ╚═╝  ╚══════╝  ╚═╝  ╚═╝  ╚═╝  ╚═════╝ \n\n"
" ---- Message from the Author! ----\n\n"
" ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓\n\n" 
"Make sure you are running the correct WORKFLOW Module.\n"
"Snakemake will automatically detect the MAIN workflow for Raw Sequence FASTQ data since the filename is -Snakefile-\n\n"
"  ►► --snakefile workflow/Snakefile ◄◄ [DEFAULT - no flag needed]\n\n"
"Snakemake requires the --snakefile flag to initiate the specific modules.\n"
"Thus, if you would like to use the SNP Array Module - Add flag below:\n\n"
"  --snakefile workflow/snpArray_Snakefile\n\n"
"Workflow Module Options:\n"
"  --snakefile workflow/cram_Snakefile\n"
"  --snakefile workflow/liftoverTo37_Snakefile\n"
"  --snakefile workflow/liftoverTo38_Snakefile\n"
"  --snakefile workflow/mergerSub_Snakefile\n"
"  --snakefile workflow/targetRefMerge_Snakefile\n\n"
"**ALL Snakefiles are all located in ./Iliad/workflow/**\n\n"
" ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑\n"
" ----------------------------------\n")

##### Helper functions #####

def get_reads(wildcards):
	"""Get condensed reads from all lanes of given sample."""
	return expand(
		"results/fastq/{{sample}}.R{group}.fq",
		group=[1, 2]
	)

def userTable_to_fqDict(fqSampleTablePath):
	#converting the .xlsx file if needed
	extension = pathlib.Path(fqSampleTablePath).suffix
	names = ["sampleID","URL"]
	fqDict = {}
	if extension == '.xlsx':
		data = pd.read_excel(fqSampleTablePath,header=None,names=names)
		data.sort_values(by=["sampleID"],axis=0, ascending=[True],inplace=True)
		data = data.dropna()
			# Creating and returning the dict 
		for row, the_value in enumerate(data['sampleID']):
			fqDict.setdefault(the_value, []).append(data.iat[row,1])
		return fqDict
	elif extension == '.csv':
		data = pd.read_csv(fqSampleTablePath,header=None,names=names)
		data.sort_values(by=["sampleID"],axis=0, ascending=[True],inplace=True)
		data = data.dropna()
			# Creating and returning the dict 
		for row, the_value in enumerate(data['sampleID']):
			fqDict.setdefault(the_value, []).append(data.iat[row,1])
		return fqDict
	return fqDict
fqDict = userTable_to_fqDict(config["samplesDict"])


def cramSampleTable_to_fqDict(cramSampleTablePath):
	#converting the .xlsx file if needed
	extension = pathlib.Path(cramSampleTablePath).suffix
	names = ["sampleID","URL"]
	cramDict = {}
	if extension == '.xlsx':
		data = pd.read_excel(cramSampleTablePath,header=None,names=names)
		data.sort_values(by=["sampleID"],axis=0, ascending=[True],inplace=True)
		data = data.dropna()
			# Creating and returning the dict 
		for row, the_value in enumerate(data['sampleID']):
			cramDict.setdefault(the_value, []).append(data.iat[row,1])
		return cramDict
	elif extension == '.csv':
		data = pd.read_csv(cramSampleTablePath,header=None,names=names)
		data.sort_values(by=["sampleID"],axis=0, ascending=[True],inplace=True)
		data = data.dropna()
			# Creating and returning the dict 
		for row, the_value in enumerate(data['sampleID']):
			cramDict.setdefault(the_value, []).append(data.iat[row,1])
		return cramDict
	return cramDict
cramDict = cramSampleTable_to_fqDict(config["cramSamplesDict"])

def which_ref(wildcards):
	input_list = []
	if config["AutoRetrieveReference"]:
		  input_list.append("resources/genome.fasta")
	if config["IhaveReference"]:
		  input_list.append(config["reference"]["filePath"])
	return input_list


# potentially provide both 37 and 38 to each rule at decision point....


def get_splits(wildcards):
	"""Get splits for all chromosomes of given sample."""
	return expand(
	"NYGC/SplitLists/chr{{chroms}}/regionsSplit-{{split}}.txt")

def get_splits_chr(wildcards):
	"""Get splits for all chromosomes of given sample."""
	return expand(
	"NYGC/SplitLists/chr{{chroms}}/regionsSplit-chr-{{split}}.txt")





def get_variant_calling_mpileup_params(wildcards):
	return config["VariantCalling"]["mpileup"]["options"]

def get_variant_calling_call_params(wildcards):
	return config["VariantCalling"]["call"]["options"]

def get_NormTrueFalse(wildcards):
	if config["Normalize"]:
		return True
	else:
		return False

def get_NormLeftAlign(wildcards):
	input_list = []
	if config["Normalize"]:
		if config["Norm"]["options"] is not None:
			input_list.append(config["Norm"]["options"])
	if config["doNotNormalize"]:
		input_list.append("")
	return input_list

def get_split_filepath(wildcards):
	"""Get split filepaths for all chromosomes of given sample."""
	return expand("results/vcf/{{sample}}/chr{{chroms}}-split-{split}-filePath.txt",
		split=[f"{i}".zfill(3) for i in range(splits)]
	)

def get_RefChr_filepath(wildcards):
	"""Get Reference Chromosome filepaths."""
	return "data/target_ref_merge/Ref-chr{wildcards.chroms}-filePath.txt"

def get_split_vcf_concat_list(wildcards):
	"""Get splits for all chromosomes of given sample."""
	split=[f"{i}".zfill(3) for i in range(splits)]
	for c in CHROMS:
		concatList=[]
		# open file in write mode
		with open(f"results/vcf/{wildcards.sample}/chr{c}-concatList.txt", 'w') as fp:
			for x in split:
				# write each vcf to concat on a new line
				addVCF=f"results/vcf/{wildcards.sample}/chr{c}.split-{x}.vcf.gz"
				fp.write("%s\n" % addVCF)
	return "results/vcf/{wildcards.sample}/chr{chroms}-concatList.txt"

def get_v37_filePath_files(wildcards):
	print(f"get_v37_filePath_files")
	print(f"get_v37_filePath_files {wildcards}")
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step2-ReferenceAssembly-VersionCheck/filePath_If_V37.{vcf}.txt", 
	vcf=vcfs["baseFileName_VCF"])

def get_v37_basenames(wildcards):
	print(f"get_v37_basenames")
	print(f"get_v37_basenames {wildcards}")
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step2-ReferenceAssembly-VersionCheck/basename_If_V37.{vcf}.txt", 
	vcf=vcfs["baseFileName_VCF"])

def get_v38_filePath_files(wildcards):
	print(f"get_v38_filePath_files")
	print(f"get_v38_filePath_files {wildcards}")
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step2-ReferenceAssembly-VersionCheck/filePath_If_V38.{vcf}.txt", 
	vcf=vcfs["baseFileName_VCF"])

def get_v38_basenames(wildcards):
	print(f"get_v38_basenames")
	print(f"get_v38_basenames {wildcards}")
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step2-ReferenceAssembly-VersionCheck/basename_If_V38.{vcf}.txt", 
	vcf=vcfs["baseFileName_VCF"])



# # input function for the rule concat_splits_per_chrom
# def get_concatenate_input(wildcards):
# 	# decision based on content of output file
# 	# Important: use the method open() of the returned file!
# 	# This way, Snakemake is able to automatically download the file if it is generated in
# 	# a cloud environment without a shared filesystem.
# 	with checkpoints.concat_list.get(sample=wildcards.sample, chroms=wildcards.chroms).output[0].open() as f:
# 		x = len(f.readlines())
# 		if x == splits:
# 			return "results/vcf/{sample}/chr{chroms}-concatList.txt"
# 		else:
# 			print(f"split vcf files have not all been added to concat list: {splits}")

# def get_versionCheck(wildcards):
# 	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step2-ReferenceAssembly-VersionCheck/versionCheck.{vcf}.txt",
# 		 vcf=vcfs["baseFileName_VCF"])
# def get_dbSNP_Extracted(wildcards):
# 	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step2-ReferenceAssembly-VersionCheck/random1.{vcf}-sites-extracted-from-dbSNP.txt",
# 		 vcf=vcfs["baseFileName_VCF"])
# def get_keep23(wildcards):
# 	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step1-Keep23/Keep23.{vcf}.vcf.gz",
# 		 vcf=vcfs["baseFileName_VCF"])
# def get_keep23index(wildcards):
# 	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step1-Keep23/Keep23.{vcf}.vcf.gz.tbi",
# 		 vcf=vcfs["baseFileName_VCF"])



# def get_combineMyData_37filepaths(wildcards):
# 	checkpoints.check_basenames.get(wildcards.project, wildcards.refAssemblyVersion)
# 	vcf37s = glob_wildcards(f"data/vcf_Merge-and-Lift/{wildcards.project}/{wildcards.refAssemblyVersion}/step4A-Input-vcfIDs37/queried.{{vcf37}}.rsIDs.txt").vcf37
# 	return expand(f"data/vcf_Merge-and-Lift/{wildcards.project}/{wildcards.refAssemblyVersion}/step4A-Input-vcfIDs37/queried.{{vcf37}}.rsIDs.txt", vcf37=vcf37s)


# def get_combineMyData_37filepaths(wildcards):
# 	"""Get combineMyData filepaths for all rsid files from merger VCFs."""
# 	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step4A-Input-vcfIDs37/queried.{vcf}.rsIDs.txt", vcf=vcfs["baseFileName_VCF"])

# def get_combineMyData_38filepaths(wildcards):
# 	"""Get combineMyData filepaths for all rsid files from merger VCFs."""
# 	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step4B-Input-vcfIDs38/queried38.{vcf}.rsIDs.txt", vcf=vcfs["baseFileName_VCF"])




# def get_reference_assembly_version(wildcards):
# 	"""Get version number for proper dbSNP liftover."""

# def get_multiQC_input(wildcards):
# 	fastQCdirs = expand("results/fastq/{sample}-fastqc/", sample=list(wildcards.sample))
# 	return fastQCdirs

