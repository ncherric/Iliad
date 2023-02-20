# Add If statement for .fastq.gz extensions or integrate into shell 

rule concat_reads:
	input:
		SMP = "results/downloads/{sample}/",
	output:
		R1="results/fastq/{sample}.R1.fq",
		R2="results/fastq/{sample}.R2.fq",
	benchmark:
		repeat("benchmarks/concat_reads.{sample}.RawSeq", 1)
	singularity:
		"docker://redhat/ubi8"
	params:
		workdir=config["workdirPath"],
	resources:
		mem_mb=25000,
		runtime="03:30:00",
	shell:
		"""
		zcat $( ls {params.workdir}{input.SMP}/*R1.fq.gz ) > {output.R1};
		zcat $( ls {params.workdir}{input.SMP}/*R2.fq.gz ) > {output.R2}
		"""