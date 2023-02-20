rule download_referenceGenome_37_file: #wont be able to fully cache as output stands
	output:
		tmpoGenomeFile37="resources/tempFile37.to.remove",
	params:
		GenomeFastaUrl=config["genomeReference"]["37Reference"],
		genomeFastaDir="resources/",
		workdirPath=config['workdirPath'],
	benchmark:
		"benchmarks/download_dbSNP_37_file.LiftAndMerge"
	resources:
		mem_mb=500,
		runtime="06:00:00",
	shell:
		"""
		mkdir -p {params.workdirPath}{params.genomeFastaDir}
		wget -t 0 --retry-connrefused -c -nH -P {params.workdirPath}{params.genomeFastaDir} {params.GenomeFastaUrl}.fasta.gz
		sleep 10s
		wget -t 0 --retry-connrefused -c -nH -P {params.workdirPath}{params.genomeFastaDir} {params.GenomeFastaUrl}.fasta.fai
		sleep 10s
		gunzip --quiet {params.workdirPath}{params.genomeFastaDir}human_g1k_v37.fasta.gz || true
		sleep 10s
		echo "to remove" > {params.workdirPath}{output.tmpoGenomeFile37}
		"""

rule download_referenceGenome_38_file: #wont be able to fully cache as output stands
	output:
		tmpoGenomeFile38="resources/tempFile38.to.remove",
	params:
		GenomeFastaUrl=config["genomeReference"]["38Reference"],
		GenomeFasta=config["genomeReference"]["file38"],
		GenomeFastaIndex=config["genomeReference"]["index38"],
		genomeFastaDir="resources/",
		workdirPath=config['workdirPath'],
	benchmark:
		"benchmarks/download_dbSNP_38_file.LiftAndMerge"
	resources:
		mem_mb=500,
		runtime="06:00:00",
	shell:
		"""
		mkdir -p {params.workdirPath}{params.genomeFastaDir}
		wget -t 0 --retry-connrefused -c -nH -P {params.workdirPath}{params.genomeFastaDir} {params.GenomeFastaUrl}{params.GenomeFasta}
		sleep 10s
		wget -t 0 --retry-connrefused -c -nH -P {params.workdirPath}{params.genomeFastaDir} {params.GenomeFastaUrl}{params.GenomeFastaIndex}
		sleep 10s
		echo "to remove" > {params.workdirPath}{output.tmpoGenomeFile38}
		"""
