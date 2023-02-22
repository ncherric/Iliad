rule download_dbSNP_file:
	output:
		tmpoFile="dbSNP/tempFile.to.remove",
	params:
		dbsnpVcfUrl=config["dbSNP"]["dbsnpVcfDownload"],
		dbsnpTbiUrl=config["dbSNP"]["dbsnpTbiDownload"],
		dbsnpDir="dbSNP/",
		workdirPath=config['workdirPath'],
	cache: True
	benchmark:
		"benchmarks/download_dbSNP_file.StoredSequence"
	resources:
		mem_mb=500,
		runtime="06:00:00",
	shell:
		"""
		mkdir -p {params.workdirPath}{params.dbsnpDir}
		wget -t 0 --retry-connrefused -c -nH -P {params.workdirPath}{params.dbsnpDir} {params.dbsnpVcfUrl}
		sleep 10s
		wget -t 0 --retry-connrefused -c -nH -P {params.workdirPath}{params.dbsnpDir} {params.dbsnpTbiUrl}
		sleep 10s
		echo "to remove" > {params.workdirPath}{output.tmpoFile}
		"""
