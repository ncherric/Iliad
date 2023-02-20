# rule fix_ref:
# 	input:
# 		filteredGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/three.Filtered.{vcf}.vcf.gz",
# 		filteredTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/three.Filtered.{vcf}.vcf.gz.tbi",
# 		# sortedGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
# 		# sortedTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz.csi",
# 		# ref=rules.bwa_index.output,
# 		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.fa",
# 		tmpoFile="dbSNP/tempFile.to.remove",
# 		# dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz",
# 		# dbsnpExtractedIDsFileTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz.tbi",
# 	output:
# 		FixRefGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step7-fixref/fixref.{vcf}.vcf.gz",
# 	params:
# 		dbsnpDir="dbSNP/",
# 		dbsnpFile=config['dbSNP']['file'],
# 		project_from_wc=lambda wc: wc.get("project"),
# 		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
# 		vcf_from_wc=lambda wc: wc.get("vcf"),
# 	resources:
# 		mem_mb=9000,
# 		runtime="01:00:00",
# 	benchmark:
# 		repeat("benchmarks/fix_ref-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	script:
# 		"../scripts/fixref.py"

# rule fix_ref_index:
# 	input:
# 		FixRefGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step7-fixref/fixref.{vcf}.vcf.gz",
# 	output:
# 		FixRefTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step7-fixref/fixref.{vcf}.vcf.gz.tbi",
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:30:00",
# 	benchmark:
# 		repeat("benchmarks/fix_ref_index-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools index --tbi {input.FixRefGZ} -o {output.FixRefTBI}
# 		"""

rule Merger_VCF_List:
	input:
		FixRefGZ=expand("data/vcf_merge/{{project}}/{{refAssemblyVersion}}/step6-liftedOver/three.Filtered.{vcf}.vcf.gz", vcf=vcfs["baseFileName_VCF"]),
		FixRefGZ_TBI=expand("data/vcf_merge/{{project}}/{{refAssemblyVersion}}/step6-liftedOver/three.Filtered.{vcf}.vcf.gz.tbi", vcf=vcfs["baseFileName_VCF"]),
		# FixRefGZ=expand("data/vcf_merge/{{project}}/{{refAssemblyVersion}}/step7-fixref/fixref.{vcf}.vcf.gz", vcf=vcfs["baseFileName_VCF"]),
		# FixRefGZ_TBI=expand("data/vcf_merge/{{project}}/{{refAssemblyVersion}}/step7-fixref/fixref.{vcf}.vcf.gz.tbi", vcf=vcfs["baseFileName_VCF"]),
	output:
		mergeList="data/vcf_merge/{project}/{refAssemblyVersion}/step8-merge/mergerList.txt",
	resources:
		mem_mb=200,
		runtime="00:02:00",
	benchmark:
		repeat("benchmarks/Merger_VCF_List-{project}-{refAssemblyVersion}.merger", 1)
	script:
		"../scripts/merger-VCF-list.py"

rule merge_vcfs:
	input:
		mergeList="data/vcf_merge/{project}/{refAssemblyVersion}/step8-merge/mergerList.txt",
	output:
		mergeVCF="data/vcf_merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz",
	params:
	resources:
		mem_mb=10000,
		runtime="00:60:00",
	benchmark:
		repeat("benchmarks/merge_vcfs-{project}-{refAssemblyVersion}.merger", 1)
	script:
		"../scripts/finalMerge.py"

rule merge_vcfs_index:
	input:
		mergeVCF="data/vcf_merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz",
	output:
		mergeTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz.tbi",
	# conda: "../../env/bcftools1-14.yaml",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/merge_vcfs_index-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools index --tbi {input.mergeVCF} -o {output.mergeTBI}
		"""
