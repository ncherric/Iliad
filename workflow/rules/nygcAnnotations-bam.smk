configfile: "config/config.yaml"

splits=config['NYGC']['numberOfSplitRegionsFiles']

checkpoint download_NYGC_annotations_file:
	output:
		AnnotationFile="NYGC/CCDG_13607_B01_GRM_WGS_2019-02-19_chr{chroms}.recalibrated_variants.annotated.txt",
	params:
		nygcUrlPath=config['NYGC']['nygcUrlPath'],
		nygcFileStart=config['NYGC']['nygcFileStart'],
		nygcFileEnd=config['NYGC']['nygcFileEnd'],
		nygcDir="NYGC/",
		splitDir="SplitLists/",
		splitLength="20000",
		numberOfSplits=config['NYGC']['numberOfSplitRegionsFiles'],
		workdirPath=config['workdirPath'],
	benchmark:
		repeat("benchmarks/download_NYGC_annotations_file.{chroms}.RawSeq", 1)
	resources:
		mem_mb=500,
		runtime="03:30:00",
	shell:
		"""
		mkdir -p {params.workdirPath}{params.nygcDir}
		wget -t 0 --retry-connrefused -c -nH -O {params.workdirPath}{output.AnnotationFile} {params.nygcUrlPath}{params.nygcFileStart}{wildcards.chroms}{params.nygcFileEnd}
		"""	

rule unique_regions_file:
	input:
		AnnotationFile="NYGC/CCDG_13607_B01_GRM_WGS_2019-02-19_chr{chroms}.recalibrated_variants.annotated.txt",
	output:
		uniq="NYGC/annotatedNYGC-chr{chroms}.uniq.txt",
	resources:
		mem_mb=10000,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/unique_regions_file.{chroms}.RawSeq", 1)
	shell:
		"""
		awk '{{ print $1,$2 }}' {input.AnnotationFile} | uniq > {output.uniq}
		"""

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

# chrStr = findSequenceName(sample=wildcards.sample)
# chrStr=lambda wildcards: findSequenceName('{sample}'.format(sample=wildcards.sample))

ruleorder: tabbed_regions_chr_string > tabbed_regions

# if chrStr:
rule tabbed_regions_chr_string: # CHROMOSOME STRING YES - "chr1"
	input:
		uniqFile="NYGC/annotatedNYGC-chr{chroms}.uniq.txt",
	output:
		tabChrFile="NYGC/annotatedNYGC-chr{chroms}.uniq.tab.chr.txt",
	resources:
		mem_mb=500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/tabbed_regions_chr_string.{chroms}.StoredSequence", 1)
	script:
		"../scripts/tabbed-regions-chrString.py"
# else:
rule tabbed_regions: # CHROMOSOME STRING NO - "1"
	input:
		uniqFile="NYGC/annotatedNYGC-chr{chroms}.uniq.txt",
	output:
		tabFile="NYGC/annotatedNYGC-chr{chroms}.uniq.tab.txt",
	resources:
		mem_mb=500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/tabbed_regions.{chroms}.StoredSequence", 1)
	script:
		"../scripts/tabbed-regions.py"

ruleorder: split_chroms_chr_string > split_chroms

rule split_chroms_chr_string:
	input:
		tabChrFile="NYGC/annotatedNYGC-chr{chroms}.uniq.tab.chr.txt",
	output:
		splitsChr="NYGC/SplitLists/chr{chroms}/regionsSplit-chr-{split}.txt",
	params:
		numberOfSplits=config['NYGC']['numberOfSplitRegionsFiles'],
		wd=config['workdirPath'],
		splitChrDir="NYGC/SplitLists/chr",
	resources:
		mem_mb=5000,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/split_chroms.chr{chroms}.split{split}.RawSeq", 1)
	shell:
		"""
		mkdir -p {params.splitChrDir}{wildcards.chroms}/
		touch {output.splitsChr}
		cd {params.splitChrDir}{wildcards.chroms}/
		split -a 3 -n l/{params.numberOfSplits} -d {params.wd}{input.tabChrFile} regionsSplit-chr- --additional-suffix=.txt
		"""

rule split_chroms:
	input:
		tabFile="NYGC/annotatedNYGC-chr{chroms}.uniq.tab.txt",
	output:
		splits="NYGC/SplitLists/chr{chroms}/regionsSplit-{split}.txt",
	params:
		numberOfSplits=config['NYGC']['numberOfSplitRegionsFiles'],
		wd=config['workdirPath'],
		splitChrDir="NYGC/SplitLists/chr",
	resources:
		mem_mb=5000,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/split_chroms.{chroms}.split{split}.RawSeq", 1)
	shell:
		"""
		mkdir -p {params.splitChrDir}{wildcards.chroms}/
		touch {output.splits}
		cd {params.splitChrDir}{wildcards.chroms}/
		split -a 3 -n l/{params.numberOfSplits} -d {params.wd}{input.tabFile} regionsSplit- --additional-suffix=.txt
		"""





