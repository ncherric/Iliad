def filtered_38_vcfs(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step7A-38_lift_to_37/QC-Filtered.{vcf38}.vcf.gz",vcf38=baseNames)

def filtered_38_index(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step7A-38_lift_to_37/QC-Filtered.{vcf38}.vcf.gz.tbi",vcf38=baseNames)


def filtered_37_vcfs(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v37_again.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step7B-no_lift_needed_37/QC-Filtered37.{valid37}.vcf.gz",valid37=baseNames)

def filtered_37_index(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.validate_basenames_v37_again.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "valid_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step7B-no_lift_needed_37/QC-Filtered37.{valid37}.vcf.gz.tbi",valid37=baseNames)


rule Merger_VCF_List:
	input:
		VCFs_in_38=filtered_38_vcfs,
		VCFs_in_37=filtered_37_vcfs,
		VCFindex_in_38=filtered_38_index,
		VCFindex_in_37=filtered_37_index,
	output:
		mergeList="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/mergeList.txt",
	params:
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
	resources:
		mem_mb=200,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/Merger_VCF_List-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-merge-VCF-list.py"

rule merge_vcfs:
	input:
		mergeList="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/mergeList.txt",
	output:
		mergeVCF="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz",
	params:
	resources:
		mem_mb=10000,
		runtime="00:60:00",
	benchmark:
		repeat("benchmarks/merge_vcfs-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-finalMerge.py"

rule merge_vcfs_index:
	input:
		mergeVCF="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz",
	output:
		mergeTBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/merge_vcfs_index-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools index --tbi {input.mergeVCF} -o {output.mergeTBI}
		"""

rule FinalQC:
	input:
		mergeVCF="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz",
		mergeTBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step8-merge/FinalMerge.vcf.gz.tbi",
	output:
		FinalQC="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step9-FinalQC/FinalMerge-QC.vcf.gz",
	params:
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/FinalQC-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i 'AC!=0 & MAF[0]<0.05 & F_MISSING < 0.05' {input.mergeVCF} -O z -o {output.FinalQC}
		"""
   
rule FinalQC_index:
	input:
		FinalQC="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step9-FinalQC/FinalMerge-QC.vcf.gz",
	output:
		FinalQC_TBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step9-FinalQC/FinalMerge-QC.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/FinalQC_index-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools index --tbi {input.FinalQC} -o {output.FinalQC_TBI}
		"""

