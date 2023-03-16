rule get_clean_rsids38_from_MyData:
	input:
		cleanDBSNPlist38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5B-dbSNP-IDs-38/queried38.clean-dbSNP-combinedMyData38.rsIDs.txt",
		annotated38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/annotated38.{valid38}.vcf.gz",
		annotated38index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4B-Input-vcfIDs38/annotated38.{valid38}.vcf.gz.tbi",
	output:
		cleanDBSNPmyData38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38.{valid38}.vcf.gz",
		cleanDBSNPmyData38Index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38.{valid38}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/get_clean_rsids38_from_MyData-{valid38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i ID=@{input.cleanDBSNPlist38} -O z -o {output.cleanDBSNPmyData38} {input.annotated38}
		bcftools index --tbi {output.cleanDBSNPmyData38}
		"""

# fixref wont work on already lifted over data just so you know! 
rule rename_chrs38:
	input:
		cleanDBSNPmyData38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38.{valid38}.vcf.gz",
		cleanDBSNPmyData38Index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38.{valid38}.vcf.gz.tbi",
	output:
		renamedChrs38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38-chrString.{valid38}.vcf.gz",
		renamedChrs38index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38-chrString.{valid38}.vcf.gz.tbi",
	params:
		renameChrFile38="config/renameChrFile38.txt",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/rename_chrs38-{valid38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate --rename-chrs {params.renameChrFile38} -O z -o {output.renamedChrs38} {input.cleanDBSNPmyData38}
		bcftools index --tbi {output.renamedChrs38}
		"""


rule fix_ref38:
	input:
		renamedChrs38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38-chrString.{valid38}.vcf.gz",
		renamedChrs38index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38-chrString.{valid38}.vcf.gz.tbi",
		tmpoGenomeFile38="resources/tempFile38.to.remove",
		tmpoFile="dbSNP/tempFile38.to.remove",
	output:
		FixRef38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/fixref38.{valid38}.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
		genomeFastaDir="resources/",
		genomeFastaFile=config['genomeReference']['file38'],
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		valid38_from_wc=lambda wc: wc.get("valid38"),
	resources:
		mem_mb=9000,
		runtime="08:00:00",
	benchmark:
		repeat("benchmarks/fix_ref38-{valid38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-fixref38.py"


# sort 
rule sort38:
	input:
		FixRef38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/fixref38.{valid38}.vcf",
	output:
		sortGZ38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38.{valid38}.vcf.gz",
		sort38TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38.{valid38}.vcf.gz.tbi",
	params:
		tempDir="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/TempDir/{valid38}/"
	resources:
		mem_mb=6000,
		runtime="01:00:00",
	benchmark:
		repeat("benchmarks/sort38-{valid38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		mkdir -p {params.tempDir}
		bcftools sort -m 5G -T {params.tempDir} -O z -o {output.sortGZ38} {input.FixRef38}
		bcftools index --tbi {output.sortGZ38}
		"""

rule rename_back_to_noChrs:
	input:
		sortGZ38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38.{valid38}.vcf.gz",
		sort38TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38.{valid38}.vcf.gz.tbi",
	output:
		noChrName38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38-noChrString.{valid38}.vcf.gz",
		noChrName38TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38-noChrString.{valid38}.vcf.gz.tbi",
	params:
		renameChrFile="config/renameChrFile37.txt",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/rename_back_to_noChrs-{valid38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate --rename-chrs {params.renameChrFile} -O z -o {output.noChrName38} {input.sortGZ38}
		bcftools index --tbi {output.noChrName38}
		"""

# qc Filtered
rule filter38:
	input:
		noChrName38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38-noChrString.{valid38}.vcf.gz",
		noChrName38TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/Sorted38-noChrString.{valid38}.vcf.gz.tbi",
	output:
		filtered38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/QC-Filtered38.{valid38}.vcf.gz",
		filtered38TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/QC-Filtered38.{valid38}.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:40:00",
	benchmark:
		repeat("benchmarks/filter38-{valid38}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools view -i 'F_MISSING < 0.05' {input.noChrName38} -O z -o {output.filtered38}
		bcftools index --tbi {output.filtered38}
		"""

# # filtering by skipping fixref step if you are attempting to use already lifted data. e.g. you already had some VCF that had positions already lifted over
# rule filter38:
# 	input:
# 		cleanDBSNPmyData38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38.{valid38}.vcf.gz",
# 		cleanDBSNPmyData38Index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/cleaned38.{valid38}.vcf.gz.tbi",
# 	output:
# 		filtered38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/QC-Filtered38.{valid38}.vcf.gz",
# 		filtered38TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7B-no_lift_needed_38/QC-Filtered38.{valid38}.vcf.gz.tbi",
# 	params:
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:40:00",
# 	benchmark:
# 		repeat("benchmarks/filter38-{valid38}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools view -i 'F_MISSING < 0.05' {input.cleanDBSNPmyData38} -O z -o {output.filtered38}
# 		bcftools index --tbi {output.filtered38}
# 		"""

def filtered_38_vcfs(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v38_again.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step7B-no_lift_needed_38/QC-Filtered38.{valid38}.vcf.gz.tbi",valid38=baseNames)

rule finished3:
	input: 
		filtered_38_vcfs,
	output:
		finished3="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/finished-step3.txt"
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	shell:
		"""
		touch {output.finished3}
		"""