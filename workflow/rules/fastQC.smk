rule FastQC_plots:
	input:
		reads=get_reads,
	output:
		directory("results/fastq/{sample}-fastqc/"),
	singularity:
		"docker://biocontainers/fastqc:v0.11.9_cv8",
	resources:
		mem_mb=5000,
		runtime="01:30:00",
	benchmark:
		repeat("benchmarks/FastQC_plots.{sample}.RawSeq", 1),
	shell:
		"""
		mkdir -p {output}
		fastqc -o {output} {input.reads}
		"""

# rule MultiQC_plots:
# 	input:
# 		fastqcResults=get_multiQC_input,
# 	output:
# 		directory("results/fastq/multiQC/"),
# 	# conda: "../../env/multiQC.yaml",
# 	singularity:
# 		"docker://ewels/multiqc",
# 	params:
# 		workdirPath=config['workdirPath'],
# 	shell:
# 		"""
# 		multiqc {params.workdirPath}{input.fastqcResults}
# 		"""