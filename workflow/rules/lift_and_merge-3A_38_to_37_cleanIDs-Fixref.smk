rule get_clean_rsids_from_MyData:
	input:
		cleanDBSNPlist="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5A-dbSNP-IDs-38/queried.clean-dbSNP-combinedMyData38.rsIDs.txt",
		annotated38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/annotated.{vcf38}.vcf.gz",
		annotated38index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/annotated.{vcf38}.vcf.gz.tbi",
	output:
		cleanDBSNPmyData="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs38/cleaned.{vcf38}.vcf.gz",
		cleanDBSNPmyDataIndex="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs38/cleaned.{vcf38}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/get_clean_rsids_from_MyData-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i ID=@{input.cleanDBSNPlist} -O z -o {output.cleanDBSNPmyData} {input.annotated38}
		bcftools index --tbi {output.cleanDBSNPmyData}
		"""

rule fix_ref:
	input:
		cleanDBSNPmyData="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs38/cleaned.{vcf38}.vcf.gz",
		cleanDBSNPmyDataIndex="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs38/cleaned.{vcf38}.vcf.gz.tbi",
		tmpoGenomeFile38="resources/tempFile38.to.remove",
		tmpoFile="dbSNP/tempFile38.to.remove",
	output:
		FixRef="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs38/fixref.{vcf38}.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
		genomeFastaDir="resources/",
		genomeFastaFile=config['genomeReference']['file38'],
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf38_from_wc=lambda wc: wc.get("vcf38"),
	resources:
		mem_mb=9000,
		runtime="08:00:00",
	benchmark:
		repeat("benchmarks/fix_ref-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-fixref.py"