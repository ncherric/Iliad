import pandas as pd
ruleorder: zip_Index_1 > zip_Index_2

rule zip_Index_1:
	input:
		unzipped="data/vcf_Lift-and-Merge/{vcf}.vcf",
	output:
		zipped="data/vcf_Lift-and-Merge/{vcf}.zipped.vcf.gz",
		zippedIndexed="data/vcf_Lift-and-Merge/{vcf}.zipped.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/zip_Index_1-{vcf}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -t ^0,chr0 -O z -o {output.zipped} {input.unzipped}
		bcftools index --tbi {output.zipped}
		"""

rule zip_Index_2:
	input:
		prezipped2="data/vcf_Lift-and-Merge/{vcf}.vcf.gz",
		prezipped2tbi="data/vcf_Lift-and-Merge/{vcf}.vcf.gz.tbi",
	output:
		zipped2="data/vcf_Lift-and-Merge/{vcf}.zipped.vcf.gz",
		zippedIndexed2="data/vcf_Lift-and-Merge/{vcf}.zipped.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/zip_Index_1-{vcf}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -t ^0,chr0 -O z -o {output.zipped2} {input.prezipped2}
		bcftools index --tbi {output.zipped2}
		"""


rule rename_chrs37:
	input:
		zipped="data/vcf_Lift-and-Merge/{vcf}.zipped.vcf.gz",
		zippedIndexed="data/vcf_Lift-and-Merge/{vcf}.zipped.vcf.gz.tbi",
	output:
		renamedChrs37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step0-RenameChrs/RenameChrs.{vcf}.vcf.gz",
		renamedChrs37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step0-RenameChrs/RenameChrs.{vcf}.vcf.gz.tbi",
	params:
		renameChrFile="config/renameChrFile37.txt",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/rename_chrs37-{vcf}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate --rename-chrs {params.renameChrFile} -O z -o {output.renamedChrs37} {input.zipped}
		bcftools index --tbi {output.renamedChrs37}
		"""

rule keep_autosomes_and_x:
	input:
		renamedChrs37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step0-RenameChrs/RenameChrs.{vcf}.vcf.gz",
		renamedChrs37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step0-RenameChrs/RenameChrs.{vcf}.vcf.gz.tbi",
	output:
		keep23="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step1-Keep23/Keep23.{vcf}.vcf.gz",
		keep23index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step1-Keep23/Keep23.{vcf}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/-{vcf}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X -O z -o {output.keep23} {input.renamedChrs37} 
		bcftools index --tbi {output.keep23}
		"""

rule get_random_vars_for_match: # still mumbo jumbo of possible 37s and 38s
	input:
		keep23="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step1-Keep23/Keep23.{vcf}.vcf.gz",
		keep23index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step1-Keep23/Keep23.{vcf}.vcf.gz.tbi",
	output:
		vcfForRandomizing="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/temp.preppedForRandom.{vcf}.vcf",
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVarsForCheck.{vcf}.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf_from_wc=lambda wc: wc.get("vcf"),

	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/check_for_version-{vcf}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-random-vars-for-match.py"


rule get_vars_from_dbSNP_for_match:
	input:
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVarsForCheck.{vcf}.txt",
		tmpoFile37="dbSNP/tempFile37.to.remove",
	output:
		dbSNP_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVars.{vcf}-sites-extracted-from-dbSNP.txt",
		dbSNP_VCF_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVars.{vcf}-sites-extracted-from-dbSNP.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/check_for_version-{vcf}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -R {input.myData_randomVarsForMatch} -O v -o {output.dbSNP_VCF_ExtractedVarsForMatch} {params.dbsnpDir}{params.dbsnpFile}
		awk '/^[^#]/ {{print $1,$2,$3}}' OFS='\t' {output.dbSNP_VCF_ExtractedVarsForMatch} > {output.dbSNP_ExtractedVarsForMatch}
		"""

rule dbSNP_match_file: 
	input:
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVarsForCheck.{vcf}.txt",
		dbSNP_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVars.{vcf}-sites-extracted-from-dbSNP.txt",
	output:
		dbSNP_Matches="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomMatchesOrBlanks.{vcf}-dbSNP-matches.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/dbSNP_match_file-{vcf}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-dbSNP-match-file-37path.py"

rule check_matches:
	input:
		myData_randomVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVarsForCheck.{vcf}.txt",
		dbSNP_ExtractedVarsForMatch="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomVars.{vcf}-sites-extracted-from-dbSNP.txt",
		dbSNP_Matches="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/randomMatchesOrBlanks.{vcf}-dbSNP-matches.txt",
		keep23="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step1-Keep23/Keep23.{vcf}.vcf.gz",
		keep23index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step1-Keep23/Keep23.{vcf}.vcf.gz.tbi",
	output:
		fileWithPathIfVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/filePath_If_V37.{vcf}.txt",
		fileWithBaseNameIfVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/basename_If_V37.{vcf}.txt",
		fileWithPathIfVersion38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/filePath_If_V38.{vcf}.txt",
		fileWithBaseNameIfVersion38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/basename_If_V38.{vcf}.txt",
	conda: "../../env/determine-vcf-version.yaml"
	params:
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf_from_wc=lambda wc: wc.get("vcf"),
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/check_matches-{vcf}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-check-matches-37path.py"

rule aggregate_v37_VCFs_basenames:
	input:
		listOf_V37_VCFs_basenames=get_v37_basenames,
	output:
		fileWithBaseNamesForAllVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/All_Version37_VCFs_basenames.txt",
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

rule aggregate_v37_VCFs:
	input:
		listOf_V37_VCFs=get_v37_filePath_files,
	output:
		fileWithPathsForAllVersion37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/All_Version37_VCFs.txt",
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


checkpoint check_basenames_v37:
	input:
		basenamesFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/All_Version37_VCFs_basenames.txt",
	output:
		tmpoFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3A-InputVCFs-37/tmpChkPtFile.tmp",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/check_basenames_v37_{project}_{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		touch {output.tmpoFile}
		"""


 ################################### ###################### DO THE SAME FOR 38

rule aggregate_v38_VCFs_basenames:
	input:
		listOf_V38_VCFs_basenames=get_v38_basenames,
	output:
		fileWithBaseNamesForAllVersion38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/All_Version38_VCFs_basenames.txt",
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

rule aggregate_v38_VCFs:
	input:
		listOf_V38_VCFs=get_v38_filePath_files,
	output:
		fileWithPathsForAllVersion38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/All_Version38_VCFs.txt",
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


checkpoint check_basenames_v38:
	input:
		basenamesFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step2-ReferenceAssembly-VersionCheck/All_Version38_VCFs_basenames.txt",
	output:
		tmpoFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step3B-InputVCFs-38/tmpChkPtFile.tmp",
	params:
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/check_basenames_v38_{project}_{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		touch {output.tmpoFile}
		"""

