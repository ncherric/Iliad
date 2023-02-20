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