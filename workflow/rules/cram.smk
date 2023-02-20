ruleorder: wget_cram > rename_cram

checkpoint wget_cram:
	output:
		smp = directory("results/downloads/{sample}"),
	benchmark:
		repeat("benchmarks/wget_cram.{sample}.StoredSequence", 1),
	resources:
		mem_mb=1000,
		runtime="03:00:00",
	params:
		cutdir=config["cramUrl"]["cutdirs"],
		urlSmp=lambda wildcards: list(cramDict[wildcards.sample]),
	shell:
		"""
		mkdir -p {output.smp};
		wget -t 0 --retry-connrefused -c -nH --cut-dirs={params.cutdir} -P {output.smp} {params.urlSmp}
		"""

rule rename_cram:
	input:
		SMP = "results/downloads/{sample}/",
	output:
		cram="results/cram/{sample}.cram",
		cramIndex="results/cram/{sample}.cram.crai",
	benchmark:
		repeat("benchmarks/rename_cram.{sample}.StoredSequence", 1),
	log:
		"logs/rename_cram.{sample}.log",
	resources:
		mem_mb=1000,
		runtime="00:05:00",
	shell:
		"""
		mkdir -p results/cram/
		for file in {input.SMP}/*
			do
			  case $file in
				*.cram) mv "$file" "{output.cram}" ;;
				*.crai) mv "$file" "{output.cramIndex}" ;;
			  esac
			done
		"""