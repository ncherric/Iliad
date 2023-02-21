# Step 2 - convert .gtc files to a VCF files
rule gtc2vcf:
	input:
		bpm="illumina/product_files/Multi-EthnicGlobal_D2.bpm",
		egt="illumina/product_files/Multi-EthnicGlobal_D1_ClusterFile.egt",
		gtc="data/snp_array/gtc/",
		fasta=which_ref,
	output:
		vcf="data/snp_array/vcf/1-vcf-from-gtc.vcf.gz",
	params:
		vcfDir="data/snp_array/vcf/",
		tempDir="data/snp_array/vcf/TempDir",
		extra_mpileup=get_variant_calling_mpileup_params,
		extra_call=get_variant_calling_call_params,
		extra_NormTrueFalse=get_NormTrueFalse,
		extra_normalize=get_NormLeftAlign,
	# benchmark:
	# 	repeat("benchmarks/gtc2vcf.array", 3)
	script:
		"../scripts/gtc2vcf.py"

rule gtc2vcf_index:
	input:
		vcf="data/snp_array/vcf/1-vcf-from-gtc.vcf.gz",
	output:
		vcfTBI="data/snp_array/vcf/1-vcf-from-gtc.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	# benchmark:
	# 	repeat("benchmarks/gtc2vcf_index.array", 3)
	shell:
		"""
		bcftools index --tbi {input.vcf} -o {output.vcfTBI}
		"""