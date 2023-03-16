rule get_clean_rsids37_from_MyData:
	input:
		cleanDBSNPlist37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-37/queried37.clean-dbSNP-combinedMyData37.rsIDs.txt",
		annotated37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/annotated37.{valid37}.vcf.gz",
		annotated37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs37/annotated37.{valid37}.vcf.gz.tbi",
	output:
		cleanDBSNPmyData37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37.{valid37}.vcf.gz",
		cleanDBSNPmyData37Index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37.{valid37}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/get_clean_rsids37_from_MyData-{valid37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i ID=@{input.cleanDBSNPlist37} -O z -o {output.cleanDBSNPmyData37} {input.annotated37}
		bcftools index --tbi {output.cleanDBSNPmyData37}
		"""

# fixref wont work on already lifted over data just so you know! 
rule rename_chrs37:
	input:
		cleanDBSNPmyData37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37.{valid37}.vcf.gz",
		cleanDBSNPmyData37Index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37.{valid37}.vcf.gz.tbi",
	output:
		renamedChrs37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37-chrString.{valid37}.vcf.gz",
		renamedChrs37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37-chrString.{valid37}.vcf.gz.tbi",
	params:
		renameChrFile37="config/renameChrFile37.txt",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/rename_chrs37-{valid37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate --rename-chrs {params.renameChrFile37} -O z -o {output.renamedChrs37} {input.cleanDBSNPmyData37}
		bcftools index --tbi {output.renamedChrs37}
		"""


rule fix_ref37:
	input:
		renamedChrs37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37-chrString.{valid37}.vcf.gz",
		renamedChrs37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37-chrString.{valid37}.vcf.gz.tbi",
		tmpoGenomeFile37="resources/tempFile37.to.remove",
		tmpoFile="dbSNP/tempFile37.to.remove",
	output:
		FixRef37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/fixref37.{valid37}.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
		genomeFastaDir="resources/",
		genomeFastaFile=config['genomeReference']['file37'],
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		valid37_from_wc=lambda wc: wc.get("valid37"),
	resources:
		mem_mb=9000,
		runtime="08:00:00",
	benchmark:
		repeat("benchmarks/fix_ref37-{valid37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-fixref37.py"


# sort 
rule sort37:
	input:
		FixRef37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/fixref37.{valid37}.vcf",
	output:
		sortGZ37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37.{valid37}.vcf.gz",
		sort37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37.{valid37}.vcf.gz.tbi",
	params:
		tempDir="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/TempDir/{valid37}/"
	resources:
		mem_mb=6000,
		runtime="01:00:00",
	benchmark:
		repeat("benchmarks/sort37-{valid37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		mkdir -p {params.tempDir}
		bcftools sort -m 5G -T {params.tempDir} -O z -o {output.sortGZ37} {input.FixRef37}
		bcftools index --tbi {output.sortGZ37}
		"""

rule rename_back_to_noChrs:
	input:
		sortGZ37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37.{valid37}.vcf.gz",
		sort37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37.{valid37}.vcf.gz.tbi",
	output:
		noChrName37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37-noChrString.{valid37}.vcf.gz",
		noChrName37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37-noChrString.{valid37}.vcf.gz.tbi",
	params:
		renameChrFile="config/renameChrFile38.txt",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/rename_back_to_noChrs-{valid37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate --rename-chrs {params.renameChrFile} -O z -o {output.noChrName37} {input.sortGZ37}
		bcftools index --tbi {output.noChrName37}
		"""

# qc Filtered
rule filter37:
	input:
		noChrName37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37-noChrString.{valid37}.vcf.gz",
		noChrName37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/Sorted37-noChrString.{valid37}.vcf.gz.tbi",
	output:
		filtered37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/QC-Filtered37.{valid37}.vcf.gz",
		filtered37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/QC-Filtered37.{valid37}.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:40:00",
	benchmark:
		repeat("benchmarks/filter37-{valid37}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools view -i 'F_MISSING < 0.05' {input.noChrName37} -O z -o {output.filtered37}
		bcftools index --tbi {output.filtered37}
		"""

# # filtering by skipping fixref step if you are attempting to use already lifted data. e.g. you already had some VCF that had positions already lifted over
# rule filter37:
# 	input:
# 		cleanDBSNPmyData37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37.{valid37}.vcf.gz",
# 		cleanDBSNPmyData37Index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs37/cleaned37.{valid37}.vcf.gz.tbi",
# 	output:
# 		filtered37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/QC-Filtered37.{valid37}.vcf.gz",
# 		filtered37TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_37/QC-Filtered37.{valid37}.vcf.gz.tbi",
# 	params:
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:40:00",
# 	benchmark:
# 		repeat("benchmarks/filter37-{valid37}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools view -i 'F_MISSING < 0.05' {input.cleanDBSNPmyData37} -O z -o {output.filtered37}
# 		bcftools index --tbi {output.filtered37}
# 		"""

def filtered_37_vcfs(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v37_again.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step7B-no_lift_needed_37/QC-Filtered37.{valid37}.vcf.gz.tbi",valid37=baseNames)

rule finished3:
	input: 
		filtered_37_vcfs,
	output:
		finished3="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/finished-step3.txt"
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	shell:
		"""
		touch {output.finished3}
		"""