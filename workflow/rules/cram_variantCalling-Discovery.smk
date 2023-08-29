CHROMS=list(range(1, 23))
CHROMS.extend(['X'])

# conditional statement in sequenceChrStringCheck.smk to either proceed with 1 or add "chr" for "chr1"
rule check_for_chr_string:
	input:
		SMP = "results/cram/{sample}.cram",
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
		"../scripts/chr-string-check-cram.py"

rule cram_variant_calling:
	input:
		cram="results/cram/{sample}.cram",
		cramIndex="results/cram/{sample}.cram.crai",
		ref = which_ref,
		chrStrCheck="results/vcf/{sample}/chrStrCheck/alignmentFileHeader.out",
		regions=CHROMS,
	output:
		splitChrVCF="results/vcf/{sample}/chr{chroms}.vcf.gz",
	params:
		sample=lambda wc: wc.get("sample"),
		chroms_from_wc=lambda wc: wc.get("chroms"),
		extra_mpileup=get_variant_calling_mpileup_params,
		extra_call=get_variant_calling_call_params,
		extra_NormTrueFalse=get_NormTrueFalse,
		extra_normalize=get_NormLeftAlign,
	resources:
		mem_mb=3500,
		runtime="24:00:00",
	benchmark:
		repeat("benchmarks/cram_variant_calling.{sample}.chr{chroms}.StoredSequence", 1)
	script:
		"../scripts/cram-variant-calling-discovery.py"


rule variant_calling_index:
	input:
		ChrVCF="results/vcf/{sample}/chr{chroms}.vcf.gz",
	output:
		ChrTBI="results/vcf/{sample}/chr{chroms}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/variant_calling_index.{sample}.chr{chroms}.StoredSequence", 1)
	shell:
		"""
		bcftools index --tbi {input.ChrVCF} -o {output.ChrTBI}
		"""


rule merge_sample_chrom_list:
	input:
		concatChrVCF=expand("results/vcf/{sample}/chr{{chroms}}.vcf.gz", sample=cramSamples["cramSample"]),
		concatChrTBI=expand("results/vcf/{sample}/chr{{chroms}}.vcf.gz.tbi", sample=cramSamples["cramSample"]),
	output:
		mergeList="results/vcf/chr{chroms}-merge-sampleList.txt",
	resources:
		mem_mb=200,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/merge_sample_chrom_list.chr{chroms}.StoredSequence", 1)
	script:
		"../scripts/merge-sample-chrom-list-discovery.py"

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
		repeat("benchmarks/merge_samples_per_chrom.chr{chroms}.StoredSequence", 1)
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
		repeat("benchmarks/merge_samples_per_chrom_index.chr{chroms}.StoredSequence", 1)
	shell:
		"""
		bcftools index --tbi {input.mergeChrVCF} -o {output.mergeChrTBI}
		"""
