rule get_38_specific_dbSNP_vcf:
	input:
		combinedSNPlist="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/combinedMyData.rsIDs.txt",
		tmpoFile38="dbSNP/tempFile38.to.remove",
	output:
		projectSpecificDBSNPvcf="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/specific-38dbSNP-for-lift.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
	resources:
		mem_mb=10000,
		disk_mb=10000,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/get_38_specific_dbSNP_vcf-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-get-38-specific-dbSNP-vcf-file-for-guidefile.py"

rule create_38_guide_file:
	input:
		projectSpecificDBSNPvcf="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/specific-38dbSNP-for-lift.vcf",
	output:
		projectSpecific38guideFile="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/38.specificguideFile-from-38dbSNP-for-lift.txt",
	params:
	resources:
		mem_mb=10000,
		disk_mb=10000,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/create_38_guide_file-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-create_38_guide_file_for37to38.py"

rule lift_37_to_38:
	input:
		projectSpecific38guideFile="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/38.specificguideFile-from-38dbSNP-for-lift.txt",
		FixRef="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/fixref.{vcf37}.vcf",
	output:
		liftedOver="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/Lifted.{vcf37}.vcf",
	params:
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
	resources:
		mem_mb=50000,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/lift_37_to_38-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-37_lift_to_38.py"


rule sort:
	input:
		liftedOver="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/Lifted.{vcf37}.vcf",
	output:
		sortGZ="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/Sorted.{vcf37}.vcf.gz",
		sortTBI="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/Sorted.{vcf37}.vcf.gz.tbi",
	params:
		tempDir="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/TempDir/{vcf37}/"
	resources:
		mem_mb=6000,
		runtime="01:00:00",
	benchmark:
		repeat("benchmarks/sort-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		mkdir -p {params.tempDir}
		bcftools sort -m 5G -T {params.tempDir} -O z -o {output.sortGZ} {input.liftedOver}
		bcftools index --tbi {output.sortGZ}
		"""

rule filter:
	input:
		sortGZ="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/Sorted.{vcf37}.vcf.gz",
		sortTBI="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/Sorted.{vcf37}.vcf.gz.tbi",
	output:
		filtered="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/QC-Filtered.{vcf37}.vcf.gz",
		filteredTBI="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step7A-37_lift_to_38/QC-Filtered.{vcf37}.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:40:00",
	benchmark:
		repeat("benchmarks/filter-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i 'F_MISSING < 0.05' {input.sortGZ} -O z -o {output.filtered}
		bcftools index --tbi {output.filtered}
		"""
def filtered_37_vcfs(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v37.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "37path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Merge-and-Lift/{{project}}/{{refAssemblyVersion}}/step7A-37_lift_to_38/QC-Filtered.{vcf37}.vcf.gz.tbi",vcf37=baseNames)

rule finished2:
	input: 
		filtered_37_vcfs
	output:
		finished2="data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/finished-step2.txt"
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	shell:
		"""
		touch {output.finished2}
		"""