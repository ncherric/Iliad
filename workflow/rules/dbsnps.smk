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
		"benchmarks/download_dbSNP_file.RawSeq"
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

rule minimize_dbsnp_file:
	input:
		tmpoFile="dbSNP/tempFile.to.remove",
	output:
		dbsnpMinFile="dbSNP/dbSNP-first-5-fields.txt",
	params:
		dbsnpFile=config['dbSNP']['file'],
		dbsnpDir="dbSNP/",
		workdirPath=config['workdirPath'],
	benchmark:
		"benchmarks/minimize_dbsnp_file.RawSeq",
	resources:
		mem_mb=5500,
		runtime="06:00:00",
	shell:
		"""
		rm -f {params.workdirPath}{input.tmpoFile}
		zcat {params.dbsnpDir}{params.dbsnpFile} | awk '{{print $1,$2,$3,$4,$5}}' OFS='\t' - > {output.dbsnpMinFile}
		"""

rule dbsnp_rsids_only:
	input:
		dbsnpMinFile="dbSNP/dbSNP-first-5-fields.txt",
	output:
		rsidList="data/snp_array/vcf/E-rsids-found-from-dbSNP-all-file.txt",
	benchmark:
		"benchmarks/dbsnp_rsids_only.RawSeq",
	resources:
		mem_mb=2500,
		runtime="03:00:00",
	cache: True
	shell:
		"""
		awk '/^[^#]/ {{print $3}}' OFS='\t' {input} > {output}
		"""
