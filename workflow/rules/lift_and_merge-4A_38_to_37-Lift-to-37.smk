rule get_38_specific_dbSNP_vcf:
	input:
		combinedSNPlist="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/combinedMyData.rsIDs.txt",
		tmpoFile38="dbSNP/tempFile37.to.remove",
	output:
		projectSpecificDBSNPvcf="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/specific-37dbSNP-for-lift.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
	resources:
		mem_mb=10000,
		disk_mb=10000,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/get_37_specific_dbSNP_vcf-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-get-37-specific-dbSNP-vcf-file-for-guidefile.py"

rule create_37_guide_file:
	input:
		projectSpecificDBSNPvcf="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/specific-37dbSNP-for-lift.vcf",
	output:
		projectSpecific37guideFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/37.specificguideFile-from-37dbSNP-for-lift.txt",
	params:
	resources:
		mem_mb=10000,
		disk_mb=10000,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/create_37_guide_file-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-create_37_guide_file_for38to37.py"

rule lift_38_to_37:
	input:
		projectSpecific37guideFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/37.specificguideFile-from-37dbSNP-for-lift.txt",
		FixRef="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs38/fixref.{vcf38}.vcf",
	output:
		liftedOver="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/Lifted.{vcf38}.vcf",
	params:
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
	resources:
		mem_mb=50000,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/lift_38_to_37-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-38_lift_to_37.py"


rule sort:
	input:
		liftedOver="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/Lifted.{vcf38}.vcf",
	output:
		sortGZ="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/Sorted.{vcf38}.vcf.gz",
		sortTBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/Sorted.{vcf38}.vcf.gz.tbi",
	params:
		tempDir="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/TempDir/{vcf38}/"
	resources:
		mem_mb=6000,
		runtime="01:00:00",
	benchmark:
		repeat("benchmarks/sort-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		mkdir -p {params.tempDir}
		bcftools sort -m 5G -T {params.tempDir} -O z -o {output.sortGZ} {input.liftedOver}
		bcftools index --tbi {output.sortGZ}
		"""

rule filter:
	input:
		sortGZ="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/Sorted.{vcf38}.vcf.gz",
		sortTBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/Sorted.{vcf38}.vcf.gz.tbi",
	output:
		filtered="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/QC-Filtered.{vcf38}.vcf.gz",
		filteredTBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step7A-38_lift_to_37/QC-Filtered.{vcf38}.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:40:00",
	benchmark:
		repeat("benchmarks/filter-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i 'F_MISSING < 0.05' {input.sortGZ} -O z -o {output.filtered}
		bcftools index --tbi {output.filtered}
		"""
def filtered_38_vcfs(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step7A-38_lift_to_37/QC-Filtered.{vcf38}.vcf.gz.tbi",vcf38=baseNames)

rule finished2:
	input: 
		filtered_38_vcfs
	output:
		finished2="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/finished-step2.txt"
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	shell:
		"""
		touch {output.finished2}
		"""