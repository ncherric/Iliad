CHROMS=list(range(1, 23))
CHROMS.extend(['X'])

# conditional statement in sequenceChrStringCheck.smk to either proceed with 1 or add "chr" for "chr1"
rule check_for_chr_string:
	input:
		sortedBam="results/sortedBam/{sample}.sorted.bam",
		bamIndex="results/sortedBam/{sample}.sorted.bam.bai",
	output:
		ChrStrCheck="results/vcf/{sample}/chrStrCheck/alignmentFileHeader.out",
	params:
		sample=lambda wc: wc.get("sample"),
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/check_for_chr_string.{sample}.RawSeq", 1)
	script:
		"../scripts/chr-string-check-bam.py"


rule variant_calling:
	input:
		sortedbam="results/sortedBam/{sample}.sorted.bam",
		bamIndex="results/sortedBam/{sample}.sorted.bam.bai",
		chrStrCheck="results/vcf/{sample}/chrStrCheck/alignmentFileHeader.out",
		ref=which_ref,
		regions=CHROMS,
	output:
		ChrVCF="results/vcf/{sample}/chr{chroms}.vcf.gz",
	resources:
		mem_mb=2500,
		runtime="24:00:00",
	params:
		sample=lambda wc: wc.get("sample"),
		extra_mpileup=get_variant_calling_mpileup_params,
		extra_call=get_variant_calling_call_params,
		extra_NormTrueFalse=get_NormTrueFalse,
		extra_normalize=get_NormLeftAlign,
	benchmark:
		repeat("benchmarks/variant_calling.{sample}.chr{chroms}.RawSeq", 1)
	script:
		"../scripts/variant-calling-discovery.py"

rule variant_calling_index:
	input:
		ChrVCF="results/vcf/{sample}/chr{chroms}.vcf.gz",
	output:
		ChrTBI="results/vcf/{sample}/chr{chroms}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/variant_calling_index.{sample}.chr{chroms}.RawSeq", 1)
	shell:
		"""
		bcftools index --tbi {input.ChrVCF} -o {output.ChrTBI}
		"""

rule merge_sample_chrom_list:
	input:
		concatChrVCF=expand("results/vcf/{sample}/chr{{chroms}}.vcf.gz", sample=samples["sample"]),
		concatChrTBI=expand("results/vcf/{sample}/chr{{chroms}}.vcf.gz.tbi", sample=samples["sample"]),
	output:
		mergeList="results/vcf/chr{chroms}-merge-sampleList.txt",
	resources:
		mem_mb=200,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/merge_sample_chrom_list.chr{chroms}.RawSeq", 1)
	script:
		"../scripts/merge-sample-chrom-list.py"

rule merge_samples_per_chrom:
	input:
		mergeList="results/vcf/chr{chroms}-merge-sampleList.txt",
	output:
		mergeChrVCF="results/vcf/Merged-chr{chroms}-allSamples.vcf.gz",
	params:
		chroms_from_wc=lambda wc: wc.get("chroms"),
	resources:
		mem_mb=10000,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/merge_samples_per_chrom.chr{chroms}.RawSeq", 1)
	script:
		"../scripts/merge-sample-chroms.py"

rule merge_samples_per_chrom_index:
	input:
		mergeChrVCF="results/vcf/Merged-chr{chroms}-allSamples.vcf.gz",
	output:
		mergeChrTBI="results/vcf/Merged-chr{chroms}-allSamples.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/merge_samples_per_chrom_index.chr{chroms}.RawSeq", 1)
	shell:
		"""
		bcftools index --tbi {input.mergeChrVCF} -o {output.mergeChrTBI}
		"""
