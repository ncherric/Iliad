# maybe put these two rules in prepare_38VCFs.smk

def aggregate_38_VCFs(wildcards):
	print(f"aggregate_38_VCFs is: {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	print(f"this is checkpoint_output: {checkpoint_output}")
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	print(f"ids is: {ids}")
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	print(f"baseNames is: {baseNames}")
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-38/{{vcf38}}.vcf.gz", vcf38=baseNames)

rule get_random_vars_for_match38: # These are all likely 38s, but just double checking...
	input:
		file38=aggregate_38_VCFs,
		# file38=lambda wildcards : aggregate_38_VCFs(wildcards),
		tmpoFile="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/tmpChkPtFile.tmp",
	output:
		vcfForRandomizing="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/temp.preppedForRandom.{vcf38}.vcf",
		myData_randomVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVarsForCheck.{vcf38}.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf38_from_wc=lambda wc: wc.get("vcf38"),
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/get_random_vars_for_match38-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-random-vars-for-match38.py"


rule get_vars_from_dbSNP_for_match38:
	input:
		myData_randomVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVarsForCheck.{vcf38}.txt",
		tmpoFile38="dbSNP/tempFile38.to.remove",
	output:
		dbSNP_ExtractedVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVars.{vcf38}-sites-extracted-from-dbSNP38.txt",
		dbSNP_VCF_ExtractedVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVars.{vcf38}-sites-extracted-from-dbSNP38.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/get_vars_from_dbSNP_for_match38-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -R {input.myData_randomVarsForMatch} -O v -o {output.dbSNP_VCF_ExtractedVarsForMatch} {params.dbsnpDir}{params.dbsnpFile}
		awk '/^[^#]/ {{print $1,$2,$3}}' OFS='\t' {output.dbSNP_VCF_ExtractedVarsForMatch} > {output.dbSNP_ExtractedVarsForMatch}
		"""

rule dbSNP_match_file38: 
	input:
		myData_randomVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVarsForCheck.{vcf38}.txt",
		dbSNP_ExtractedVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVars.{vcf38}-sites-extracted-from-dbSNP38.txt",
	output:
		dbSNP_Matches="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomMatchesOrBlanks.{vcf38}-dbSNP-matches.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/dbSNP_match_file38-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-dbSNP-match-file-38path.py"

rule validate_check_matches:
	input:
		myData_randomVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVarsForCheck.{vcf38}.txt",
		dbSNP_ExtractedVarsForMatch="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomVars.{vcf38}-sites-extracted-from-dbSNP38.txt",
		dbSNP_Matches="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/randomMatchesOrBlanks.{vcf38}-dbSNP-matches.txt",
		file38=aggregate_38_VCFs,
	output:
		fileWithPathIfVersionNotValidated="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/filePath_If_not_validated-for-v37-or-v38.{vcf38}.txt",
		fileWithBaseNameIfVersionNotValidated="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/basename_If_not_validated-for-v37-or-v38.{vcf38}.txt",
		fileWithPathIfVersion38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/filePath_If_V38.{vcf38}.txt",
		fileWithBaseNameIfVersion38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/basename_If_V38.{vcf38}.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf38_from_wc=lambda wc: wc.get("vcf38"),
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/validate_check_matches-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-check-matches-38path.py"


def get_v38_filePath_files_again(wildcards):
	print(f"get_v38_filePath_files_validated")
	print(f"get_v38_filePath_files_validated {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-38/ValidateVersion/filePath_If_V38.{vcf38}.txt", vcf38=baseNames) 

def get_v38_basenames_again(wildcards):
	print(f"get_v38_basenames_validated")
	print(f"get_v38_basenames_validated {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-38/ValidateVersion/basename_If_V38.{vcf38}.txt", vcf38=baseNames)


rule aggregate_v38_VCFs_basenames_again:
	input:
		listOf_V38_VCFs_basenames=get_v38_basenames_again,
	output:
		fileWithBaseNamesForAllVersion38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/All_Version38_VCFs_basenames.txt",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/aggregate_v38_VCFs_basenames_{project}_{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input.listOf_V38_VCFs_basenames} > {output.fileWithBaseNamesForAllVersion38}
		"""

rule aggregate_v38_VCFs_again: # make this a rule? # THIS ONE GETS WHOLE FILEPATHS
	input:
		listOf_V38_VCFs=get_v38_filePath_files_again,
	output:
		fileWithPathsForAllVersion38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/All_Version38_VCFs.txt",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/aggregate_v38_VCFs-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input.listOf_V38_VCFs} > {output.fileWithPathsForAllVersion38}
		"""


checkpoint validate_basenames_v38_again:
	input:
		basenamesFile="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/ValidateVersion/All_Version38_VCFs_basenames.txt",
	output:
		tmpoFile="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/validated/tmpChkPtFile.tmp",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/validate_basenames_v38_again-{project}_{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		touch {output.tmpoFile}
		"""

def aggregate_38_VCFs_again(wildcards):
	print(f"aggregate_38_VCFs is: {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v38_again.get(**wildcards).output[0])
	print(f"this is checkpoint_output: {checkpoint_output}")
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	print(f"ids is: {ids}")
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	print(f"baseNames is: {baseNames}")
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-38/validated/{{valid38}}.vcf.gz", valid38=baseNames)



# maybe put below in this current file

rule annotate_if_VCF_38version:
	input:
		valid38=aggregate_38_VCFs_again,
		# valid38=lambda wildcards : aggregate_38_VCFs_again(wildcards),
		tmpoFile38="dbSNP/tempFile38.to.remove",
		tmpoFile="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/validated/tmpChkPtFile.tmp",
	output:
		annotated38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/annotated38.{valid38}.vcf.gz",
		annotated38index="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/annotated38.{valid38}.vcf.gz.tbi",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
	resources:
		mem_mb=3500,
		runtime="04:30:00",
	benchmark:
		repeat("benchmarks/annotate_if_VCF_38version-{valid38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate -a {params.dbsnpDir}{params.dbsnpFile} -c ID -O z -o {output.annotated38} {input.valid38[0]}
		bcftools index --tbi {output.annotated38}
		"""

rule query_38_VCF_IDs:
	input:
		annotated38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/annotated38.{valid38}.vcf.gz",
		annotated38index="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/annotated38.{valid38}.vcf.gz.tbi",
	output:
		rsidList38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/queried38.{valid38}.rsIDs.txt"
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query_38_VCF_IDs-{valid38}-{{project}}-{{refAssemblyVersion}}.LiftAndMerge", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.annotated38} > {output.rsidList38}
		"""


def get_combineMyData_38filepaths(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v38_again.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step4B-Input-vcfIDs38/queried38.{valid38}.rsIDs.txt",valid38=baseNames)

rule combine_MyData_SNPs_38:
	input:
		rsidList=get_combineMyData_38filepaths,
	output:
		combinedSNPlist38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/combinedMyData.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/combine_MyData_SNPs_38-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input} | sort -V - | uniq - > {output.combinedSNPlist38}
		"""

# maybe put below in 38dbSNP-IDs.smk

rule get_rsids_from_dbSNP38:
	input:
		combinedSNPlist38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/combinedMyData.rsIDs.txt",
	output:
		dbsnpExtractedIDsFile38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-38/extracted38.dbSNP.combinedMyData.vcf.gz",
		dbsnpExtractedIDsFile38TBI="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-38/extracted38.dbSNP.combinedMyData.vcf.gz.tbi",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/get_rsids_from_dbSNP38-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i ID=@{input.combinedSNPlist38} -O z -o {output.dbsnpExtractedIDsFile38} {params.dbsnpDir}{params.dbsnpFile}
		bcftools index --tbi {output.dbsnpExtractedIDsFile38}
		"""

rule query_dbSNP_clean_IDs38:
	input:
		dbsnpExtractedIDsFile38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-38/extracted38.dbSNP.combinedMyData.vcf.gz",
		dbsnpExtractedIDsFile38TBI="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-38/extracted38.dbSNP.combinedMyData.vcf.gz.tbi",
	output:
		cleanDBSNPlist38="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-38/queried38.clean-dbSNP-combinedMyData38.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query_dbSNP_clean_IDs38-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.dbsnpExtractedIDsFile38} > {output.cleanDBSNPlist38}
		"""
