rule annotate_rsIDs_with_dbSNP:
	input:
		mergeChrVCF="results/vcf/Merged-chr{chroms}-allSamples.vcf.gz",
		mergeChrTBI="results/vcf/Merged-chr{chroms}-allSamples.vcf.gz.tbi",
		tmpoFile="dbSNP/tempFile.to.remove",
	output:
		mergeChrVCF_rsID="results/vcf/Merged-chr{chroms}-allSamples-with-rsIDs.vcf.gz",
	params:
		dbsnpFile=config['dbSNP']['file'],
		dbsnpDir="dbSNP/",
		workdirPath=config['workdirPath'],
	resources:
		mem_mb=15000,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/annotate_rsIDs_with_dbSNP.chr{chroms}.RawSeq", 1)
	script:
		"../scripts/annotate-rsIDs-for-seq-data-vcfs.py"

rule annotation_index:
	input:
		mergeChrVCF_rsID="results/vcf/Merged-chr{chroms}-allSamples-with-rsIDs.vcf.gz",
	output:
		mergeChrVCF_rsID_TBI="results/vcf/Merged-chr{chroms}-allSamples-with-rsIDs.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/annotation_index.chr{chroms}.RawSeq", 1)
	shell:
		"""
		bcftools index --tbi {input.mergeChrVCF_rsID} -o {output.mergeChrVCF_rsID_TBI}
		"""
