rule get_rsids_from_dbSNP:
	input:
		combinedSNPlist="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/combinedMyData.rsIDs.txt",
	output:
		dbsnpExtractedIDsFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5A-dbSNP-IDs-38/extracted.dbSNP.combinedMyData.vcf.gz",
		dbsnpExtractedIDsFileTBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5A-dbSNP-IDs-38/extracted.dbSNP.combinedMyData.vcf.gz.tbi",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/get_rsids_from_dbSNP-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools view -i ID=@{input.combinedSNPlist} -O z -o {output.dbsnpExtractedIDsFile} {params.dbsnpDir}{params.dbsnpFile}
		bcftools index --tbi {output.dbsnpExtractedIDsFile}
		"""

rule query_dbSNP_clean_IDs:
	input:
		dbsnpExtractedIDsFile="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5A-dbSNP-IDs-38/extracted.dbSNP.combinedMyData.vcf.gz",
		dbsnpExtractedIDsFileTBI="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5A-dbSNP-IDs-38/extracted.dbSNP.combinedMyData.vcf.gz.tbi",
	output:
		cleanDBSNPlist="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step5A-dbSNP-IDs-38/queried.clean-dbSNP-combinedMyData38.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query_dbSNP_clean_IDs-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.dbsnpExtractedIDsFile} > {output.cleanDBSNPlist}
		"""



