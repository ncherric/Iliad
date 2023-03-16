def aggregate_37_VCFs(wildcards):
	print(f"aggregate_37_VCFs is: {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v37.get(**wildcards).output[0])
	print(f"this is checkpoint_output: {checkpoint_output}")
	ids = glob(join(checkpoint_output, "37path_*.vcf.gz"))
	print(f"ids is: {ids}")
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	print(f"baseNames is: {baseNames}")
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-37/{{vcf37}}.vcf.gz", vcf37=baseNames)

rule get_random_vars_for_match37: # These are all likely 37s, but just double checking...
	input:
		file37=aggregate_37_VCFs,
		tmpoFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/tmpChkPtFile.tmp",
	output:
		vcfForRandomizing="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/temp.preppedForRandom.{vcf37}.vcf",
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVarsForCheck.{vcf37}.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf37_from_wc=lambda wc: wc.get("vcf37"),
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/get_random_vars_for_match37-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-random-vars-for-match37.py"


rule get_vars_from_dbSNP_for_match37:
	input:
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVarsForCheck.{vcf37}.txt",
		tmpoFile37="dbSNP/tempFile37.to.remove",
	output:
		dbSNP_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVars.{vcf37}-sites-extracted-from-dbSNP37.txt",
		dbSNP_VCF_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVars.{vcf37}-sites-extracted-from-dbSNP37.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/get_vars_from_dbSNP_for_match37-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -R {input.myData_randomVarsForMatch} -O v -o {output.dbSNP_VCF_ExtractedVarsForMatch} {params.dbsnpDir}{params.dbsnpFile}
		awk '/^[^#]/ {{print $1,$2,$3}}' OFS='\t' {output.dbSNP_VCF_ExtractedVarsForMatch} > {output.dbSNP_ExtractedVarsForMatch}
		"""

rule dbSNP_match_file37: 
	input:
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVarsForCheck.{vcf37}.txt",
		dbSNP_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVars.{vcf37}-sites-extracted-from-dbSNP37.txt",
	output:
		dbSNP_Matches="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomMatchesOrBlanks.{vcf37}-dbSNP-matches.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/dbSNP_match_file37-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-dbSNP-match-file-37path.py"

rule validate_check_matches:
	input:
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVarsForCheck.{vcf37}.txt",
		dbSNP_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomVars.{vcf37}-sites-extracted-from-dbSNP37.txt",
		dbSNP_Matches="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/randomMatchesOrBlanks.{vcf37}-dbSNP-matches.txt",
		file37=aggregate_37_VCFs,
	output:
		fileWithPathIfVersionNotValidated="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/filePath_If_not_validated-for-v38-or-v37.{vcf37}.txt",
		fileWithBaseNameIfVersionNotValidated="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/basename_If_not_validated-for-v38-or-v37.{vcf37}.txt",
		fileWithPathIfVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/filePath_If_V37.{vcf37}.txt",
		fileWithBaseNameIfVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/basename_If_V37.{vcf37}.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf37_from_wc=lambda wc: wc.get("vcf37"),
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/validate_check_matches-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-check-matches-37path.py"


def get_v37_filePath_files_again(wildcards):
	print(f"get_v37_filePath_files_validated")
	print(f"get_v37_filePath_files_validated {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v37.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "37path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-37/ValidateVersion/filePath_If_V37.{vcf37}.txt", vcf37=baseNames) 

def get_v37_basenames_again(wildcards):
	print(f"get_v37_basenames_validated")
	print(f"get_v37_basenames_validated {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v37.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "37path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-37/ValidateVersion/basename_If_V37.{vcf37}.txt", vcf37=baseNames)


rule aggregate_v37_VCFs_basenames_again:
	input:
		listOf_V37_VCFs_basenames=get_v37_basenames_again,
	output:
		fileWithBaseNamesForAllVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/All_Version37_VCFs_basenames.txt",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/aggregate_v37_VCFs_basenames_{project}_{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input.listOf_V37_VCFs_basenames} > {output.fileWithBaseNamesForAllVersion37}
		"""

rule aggregate_v37_VCFs_again:
	input:
		listOf_V37_VCFs=get_v37_filePath_files_again,
	output:
		fileWithPathsForAllVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/All_Version37_VCFs.txt",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/aggregate_v37_VCFs-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input.listOf_V37_VCFs} > {output.fileWithPathsForAllVersion37}
		"""


checkpoint validate_basenames_v37_again:
	input:
		basenamesFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/ValidateVersion/All_Version37_VCFs_basenames.txt",
	output:
		tmpoFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/validated/tmpChkPtFile.tmp",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/validate_basenames_v37_again-{project}_{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		touch {output.tmpoFile}
		"""

def aggregate_37_VCFs_again(wildcards):
	print(f"aggregate_37_VCFs is: {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v37_again.get(**wildcards).output[0])
	print(f"this is checkpoint_output: {checkpoint_output}")
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	print(f"ids is: {ids}")
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	print(f"baseNames is: {baseNames}")
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step3B-InputVCFs-37/validated/{{valid37}}.vcf.gz", valid37=baseNames)

rule annotate_if_VCF_37version:
	input:
		valid37=aggregate_37_VCFs_again,
		tmpoFile37="dbSNP/tempFile37.to.remove",
		tmpoFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-37/validated/tmpChkPtFile.tmp",
	output:
		annotated37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/annotated37.{valid37}.vcf.gz",
		annotated37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/annotated37.{valid37}.vcf.gz.tbi",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
	resources:
		mem_mb=3500,
		runtime="04:30:00",
	benchmark:
		repeat("benchmarks/annotate_if_VCF_37version-{valid37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate -a {params.dbsnpDir}{params.dbsnpFile} -c ID -O z -o {output.annotated37} {input.valid37[0]}
		bcftools index --tbi {output.annotated37}
		"""

rule query_37_VCF_IDs:
	input:
		annotated37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/annotated37.{valid37}.vcf.gz",
		annotated37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/annotated37.{valid37}.vcf.gz.tbi",
	output:
		rsidList37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/queried37.{valid37}.rsIDs.txt"
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query_37_VCF_IDs-{valid37}-{{project}}-{{refAssemblyVersion}}.LiftAndMerge", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.annotated37} > {output.rsidList37}
		"""

def get_combineMyData_37filepaths(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v37_again.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step4B-Input-vcfIDs37/queried37.{valid37}.rsIDs.txt",valid37=baseNames)

rule combine_MyData_SNPs_37:
	input:
		rsidList=get_combineMyData_37filepaths,
	output:
		combinedSNPlist37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/combinedMyData.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/combine_MyData_SNPs_37-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input} | sort -V - | uniq - > {output.combinedSNPlist37}
		"""

rule get_rsids_from_dbSNP37:
	input:
		combinedSNPlist37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/combinedMyData.rsIDs.txt",
	output:
		dbsnpExtractedIDsFile37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-37/extracted37.dbSNP.combinedMyData.vcf.gz",
		dbsnpExtractedIDsFile37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-37/extracted37.dbSNP.combinedMyData.vcf.gz.tbi",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/get_rsids_from_dbSNP37-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i ID=@{input.combinedSNPlist37} -O z -o {output.dbsnpExtractedIDsFile37} {params.dbsnpDir}{params.dbsnpFile}
		bcftools index --tbi {output.dbsnpExtractedIDsFile37}
		"""

rule query_dbSNP_clean_IDs37:
	input:
		dbsnpExtractedIDsFile37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-37/extracted37.dbSNP.combinedMyData.vcf.gz",
		dbsnpExtractedIDsFile37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-37/extracted37.dbSNP.combinedMyData.vcf.gz.tbi",
	output:
		cleanDBSNPlist37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-37/queried37.clean-dbSNP-combinedMyData37.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query_dbSNP_clean_IDs37-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.dbsnpExtractedIDsFile37} > {output.cleanDBSNPlist37}
		"""
