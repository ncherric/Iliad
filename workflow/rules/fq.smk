SRA=config["SRA"]

if SRA is not None:
	print("SRA was True, retrieving fastq data from Sequence Read Archive")

	checkpoint wget_reads_from_SRA:
		output:
			smp = temp(directory("results/downloads/{sample}")),
		benchmark:
			repeat("benchmarks/wget_reads_from_SRA.{sample}.RawSeq", 1),
		params:
			# cutdir=config["url"]["cutdirs"],
			# urlSmp=lambda wildcards: list(fqDict[wildcards.sample])
		resources:
			mem_mb=1000,
			runtime="03:00:00",
		conda: "../../env/sratoolkit.yaml",
		threads: 8
		shell:
			"""
			mkdir -p {output.smp};
			fasterq-dump {wildcards.sample} -O {output.smp} -t {output.smp} -e {threads}
			"""


if SRA is None:
	print("SRA was False, retrieving fastq data from FTP site denoted in UserSampleTable.xlsx")

	checkpoint wget_reads_from_url:
		output:
			smp = temp(directory("results/downloads/{sample}")),
		benchmark:
			repeat("benchmarks/wget_reads_from_url.{sample}.RawSeq", 1),
		params:
			cutdir=config["url"]["cutdirs"],
			urlSmp=lambda wildcards: list(fqDict[wildcards.sample])
		resources:
			mem_mb=1000,
			runtime="03:00:00",
		shell:
			"""
			mkdir -p {output.smp};
			wget -t 0 --retry-connrefused -c -nH --cut-dirs={params.cutdir} -P {output.smp} {params.urlSmp}
			"""