rule download_dbSNP_37_file: #wont be able to fully cache as output stands
	output:
		tmpoFile37="dbSNP/tempFile37.to.remove",
	params:
		dbsnpVcfUrl=config["dbsnpLiftMerge"]["dbsnp37VcfDownload"],
		dbsnpTbiUrl=config["dbsnpLiftMerge"]["dbsnp37TbiDownload"],
		dbsnpDir="dbSNP/",
		workdirPath=config['workdirPath'],
	benchmark:
		"benchmarks/download_dbSNP_37_file.LiftAndMerge"
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
		echo "to remove" > {params.workdirPath}{output.tmpoFile37}
		"""

rule download_dbSNP_38_file: #wont be able to fully cache as output stands
	output:
		tmpoFile38="dbSNP/tempFile38.to.remove",
	params:
		dbsnpVcfUrl=config["dbsnpLiftMerge"]["dbsnp38VcfDownload"],
		dbsnpTbiUrl=config["dbsnpLiftMerge"]["dbsnp38TbiDownload"],
		dbsnpDir="dbSNP/",
		workdirPath=config['workdirPath'],
	benchmark:
		"benchmarks/download_dbSNP_38_file.LiftAndMerge"
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
		echo "to remove" > {params.workdirPath}{output.tmpoFile38}
		"""
