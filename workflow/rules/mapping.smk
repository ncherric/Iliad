rule alignment:
	input:
		reads=get_reads,
		ref=rules.bwa_index.output,
	output:
		sortedBam="results/sortedBam/{sample}.sorted.bam",
		bamIndex="results/sortedBam/{sample}.sorted.bam.bai",
	params:
		index=lambda w, input: os.path.splitext(input.ref[0])[0],
		sort="picard",
		sort_order="coordinate",
	benchmark:
		repeat("benchmarks/alignment.{sample}.RawSeq", 1)
	threads: 12
	resources:
		mem_mb=50000,
		runtime="72:00:00",
	# script:
	# 	"../scripts/mapping.py" # This is what works for " bwa mem -> picard tools "
	script:
		"../scripts/mapping-samtools-sort-test.py"