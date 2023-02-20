rule extract_snpArray_rsID_from_refData:
	input:
		RefChrVCF="results/vcf/Merged-chr{chroms}-allSamples.vcf.gz",
		RefChrTBI="results/vcf/Merged-chr{chroms}-allSamples.vcf.gz.tbi",
		keeperSNPs="data/snp_array/qc/SNPs-above-GenTrain-and-ClusterSep-Threshold-rsidLIST.txt",
	output:
		RefChrVCFoverlap="data/target_ref_merge/snpArray-overlap-in-RefData-chr{chroms}-allSamples.vcf.gz",
	resources:
		mem_mb=1500,
		runtime="00:60:00",
	benchmark:
		repeat("benchmarks/extract_snpArray_rsID_from_refData.chr{chroms}.mergeTargetRef", 1)
	script:
		"../scripts/extract-snpArray-rsID-from-refData.py"

rule index_snpArray_rsID_from_refData:
	input:
		RefChrVCFoverlap="data/target_ref_merge/snpArray-overlap-in-RefData-chr{chroms}-allSamples.vcf.gz",
	output:
		RefChrVCFoverlapTBI="data/target_ref_merge/snpArray-overlap-in-RefData-chr{chroms}-allSamples.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/index_snpArray_rsID_from_refData.chr{chroms}.mergeTargetRef", 1)
	shell:
		"""
		bcftools index --tbi {input.RefChrVCFoverlap} -o {output.RefChrVCFoverlapTBI}
		"""

rule ref_chr_fnPath_to_file:
	input:
		RefChrVCFoverlap="data/target_ref_merge/snpArray-overlap-in-RefData-chr{chroms}-allSamples.vcf.gz",
		RefChrVCFoverlapTBI="data/target_ref_merge/snpArray-overlap-in-RefData-chr{chroms}-allSamples.vcf.gz.tbi",
	output:
		fnPathFile="data/target_ref_merge/Ref-chr{chroms}-filePath.txt",
	resources:
		mem_mb=200,
		runtime="00:02:00",
	benchmark:
		repeat("benchmarks/ref_chr_fnPath_to_file.chr{chroms}.mergeTargetRef", 1)
	shell:
		"""
		echo "{input.RefChrVCFoverlap}" > {output}
		"""

rule cat_RefChr_list:
	input:
		RefChrFilePath=get_RefChr_filepath,
	output:
		concatList="data/target_ref_merge/Ref-concatList.txt",
	resources:
		mem_mb=200,
		runtime="00:02:00",
	benchmark:
		repeat("benchmarks/cat_RefChr_list.mergeTargetRef", 1)
	shell:
		"""
		cat {input} > {output}
		"""

rule concat_RefChr:
	input:
		concatList="data/target_ref_merge/Ref-concatList.txt",
	output:
		concatVCF="data/target_ref_merge/concatenated-Reference.vcf.gz",
	resources:
		mem_mb=10000,
		runtime="00:20:00",
	benchmark:
		repeat("benchmarks/concat_RefChr.mergeTargetRef", 1)
	script:
		"../scripts/concatenate-reference-set-overlap.py"

rule concat_RefChr_index:
	input:
		concatVCF="data/target_ref_merge/concatenated-Reference.vcf.gz",
	output:
		concatTBI="data/target_ref_merge/concatenated-Reference.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/concat_RefChr_index.mergeTargetRef", 1)
	shell:
		"""
		bcftools index --tbi {input.concatVCF} -o {output.concatTBI}
		"""

rule merge_target_and_ref:
	input:
		ReferenceVCF="data/target_ref_merge/concatenated-Reference.vcf.gz",
		TargetVCF="data/snp_array/vcf/5-passing-QC-rsids.vcf.gz",
	output:
		MergedTargetRef="data/target_ref_merge/MERGED-TARGET-REFERENCE.vcf.gz",
	resources:
		mem_mb=4500,
		runtime="03:00:00",
	benchmark:
		repeat("benchmarks/merge_target_and_ref.mergeTargetRef", 1)
	script:
		"../scripts/merge-ref-and-target.py"	

rule merged_target_reference_index:
	input:
		MergedTargetRef="data/target_ref_merge/MERGED-TARGET-REFERENCE.vcf.gz",
	output:
		MergedTargetRefTBI="data/target_ref_merge/MERGED-TARGET-REFERENCE.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/merged_target_reference_index.mergeTargetRef", 1)
	shell:
		"""
		bcftools index --tbi {input.MergedTargetRef} -o {output.MergedTargetRefTBI}
		"""
