rule get_genome:
	output:
		"resources/genome.fasta",
	log:
		"logs/get-genome.log",
	params:
		species=config["ref"]["species"],
		datatype="dna",
		build=config["ref"]["build"],
		release=config["ref"]["release"],
	resources:
		mem_mb=36900,
		runtime="06:00:00",
	benchmark:
		"benchmarks/get_genome.RawSeq",
	cache: True
	script:
		"../scripts/reference-genome.py"


rule bwa_index:
	input:
		"resources/genome.fasta",
	output:
		multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
	log:
		"logs/bwa_index.log",
	resources:
		mem_mb=36900,
		runtime="06:00:00",
	benchmark:
		"benchmarks/bwa_index.RawSeq",
	cache: True
	script:
		"../scripts/index-reference-genome-bwa.py"