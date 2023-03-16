rule get_clean_rsids_from_MyData:
	input:
		cleanDBSNPlist="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5A-dbSNP-IDs-37/queried.clean-dbSNP-combinedMyData37.rsIDs.txt",
		annotated37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/annotated.{vcf37}.vcf.gz",
		annotated37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/annotated.{vcf37}.vcf.gz.tbi",
	output:
		cleanDBSNPmyData="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/cleaned.{vcf37}.vcf.gz",
		cleanDBSNPmyDataIndex="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/cleaned.{vcf37}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/get_clean_rsids_from_MyData-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i ID=@{input.cleanDBSNPlist} -O z -o {output.cleanDBSNPmyData} {input.annotated37}
		bcftools index --tbi {output.cleanDBSNPmyData}
		"""

rule fix_ref:
	input:
		cleanDBSNPmyData="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/cleaned.{vcf37}.vcf.gz",
		cleanDBSNPmyDataIndex="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/cleaned.{vcf37}.vcf.gz.tbi",
		tmpoGenomeFile37="resources/tempFile37.to.remove",
		tmpoFile="dbSNP/tempFile37.to.remove",
	output:
		FixRef="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/fixref.{vcf37}.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
		genomeFastaDir="resources/",
		genomeFastaFile=config['genomeReference']['file37'],
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf37_from_wc=lambda wc: wc.get("vcf37"),
	resources:
		mem_mb=9000,
		runtime="08:00:00",
	benchmark:
		repeat("benchmarks/fix_ref-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	script:
		"../scripts/Lift-and-Merge-fixref.py"