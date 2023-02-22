splits=config['NYGC']['numberOfSplitRegionsFiles']

rule cram_variant_calling:
	input:
		cram="results/cram/{sample}.cram",
		cramIndex="results/cram/{sample}.cram.crai",
		ref = which_ref,
		chrStrCheck="results/vcf/{sample}/chrStrCheck/alignmentFileHeader.out",
		regionsFileChr=get_splits_chr,
		regionsFile=get_splits,
	output:
		splitChrVCF="results/vcf/{sample}/chr{chroms}.split-{split}.vcf.gz",
	params:
		sample=lambda wc: wc.get("sample"),
		chroms_from_wc=lambda wc: wc.get("chroms"),
		split_from_wc=lambda wc: wc.get("split"),
		extra_mpileup=get_variant_calling_mpileup_params,
		extra_call=get_variant_calling_call_params,
		extra_NormTrueFalse=get_NormTrueFalse,
		extra_normalize=get_NormLeftAlign,
	resources:
		mem_mb=3500,
		runtime="24:00:00",
	benchmark:
		repeat("benchmarks/cram_variant_calling.{sample}.chr{chroms}.split{split}.StoredSequence", 1)
	script:
		"../scripts/cram-variant-calling-conditional.py"


rule variant_calling_index:
	input:
		splitChrVCF="results/vcf/{sample}/chr{chroms}.split-{split}.vcf.gz",
	output:
		splitChrTBI="results/vcf/{sample}/chr{chroms}.split-{split}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/variant_calling_index.{sample}.chr{chroms}.split{split}.StoredSequence", 1)
	shell:
		"""
		bcftools index --tbi {input.splitChrVCF} -o {output.splitChrTBI}
		"""

rule fnPath_to_file:
	input:
		splitChrVCF="results/vcf/{sample}/chr{chroms}.split-{split}.vcf.gz",
		splitChrTBI="results/vcf/{sample}/chr{chroms}.split-{split}.vcf.gz.tbi",
	output:
		fnPathFile="results/vcf/{sample}/chr{chroms}-split-{split}-filePath.txt",
	resources:
		mem_mb=200,
		runtime="00:02:00",
	benchmark:
		repeat("benchmarks/fnPath_to_file.{sample}.chr{chroms}.split{split}.StoredSequence", 1)
	shell:
		"""
		echo "{input.splitChrVCF}" > {output}
		"""

rule cat_splits_to_chr_list:
	input:
		splitFilePath=get_split_filepath,
	output:
		concatList="results/vcf/{sample}/chr{chroms}-concatList.txt",
	resources:
		mem_mb=200,
		runtime="00:02:00",
	benchmark:
		repeat("benchmarks/cat_splits_to_chr_list.{sample}.chr{chroms}.StoredSequence", 1)
	shell:
		"""
		cat {input} > {output}
		"""

rule concat_splits_per_chrom:
	input:
		concatList="results/vcf/{sample}/chr{chroms}-concatList.txt",
	output:
		concatChrVCF="results/vcf/{sample}/concatenated-chr{chroms}.vcf.gz",
	params:
		samples_from_wc=lambda wc: wc.get("sample"),
		chroms_from_wc=lambda wc: wc.get("chroms"),
	resources:
		mem_mb=10000,
		runtime="00:20:00",
	benchmark:
		repeat("benchmarks/concat_splits_per_chrom.{sample}.chr{chroms}.StoredSequence", 1)
	script:
		"../scripts/concatenate-sample-chrom-splits.py"

rule concat_splits_per_chrom_index:
	input:
		concatChrVCF="results/vcf/{sample}/concatenated-chr{chroms}.vcf.gz",
	output:
		concatChrTBI="results/vcf/{sample}/concatenated-chr{chroms}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/concat_splits_per_chrom_index.{sample}.chr{chroms}.StoredSequence", 1)
	shell:
		"""
		bcftools index --tbi {input.concatChrVCF} -o {output.concatChrTBI}
		"""

rule merge_sample_chrom_list:
	input:
		concatChrVCF=expand("results/vcf/{sample}/concatenated-chr{{chroms}}.vcf.gz", sample=cramSamples["cramSample"]),
		concatChrTBI=expand("results/vcf/{sample}/concatenated-chr{{chroms}}.vcf.gz.tbi", sample=cramSamples["cramSample"]),
	output:
		mergeList="results/vcf/chr{chroms}-merge-sampleList.txt",
	resources:
		mem_mb=200,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/merge_sample_chrom_list.chr{chroms}.StoredSequence", 1)
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
