checkpoint download_MEGA_physical_genetic_coordinates:
	output:
		physicalCoordinateFile=protected("illumina/support_files/Multi-EthnicGlobal_D2.csv_Physical-and-Genetic-Coordinates.txt"),
	params:
		supportDir="illumina/support_files/",
		coord=config["urlSupportFiles"]["physicalGeneticCoordinates"],
		pzip=config["urlSupportFiles"]["pzip"],
	benchmark:
		repeat("benchmarks/download_MEGA_physical_genetic_coordinates.array", 1)
	shell:
		"""
		mkdir -p {params.supportDir}
		wget -t 0 --retry-connrefused -c -nH -P {params.supportDir} {params.coord}
		sleep 30s
		unzip -o {params.supportDir}{params.pzip} -d {params.supportDir}
		sleep 30s
		rm {params.supportDir}{params.pzip}
		"""

rule cleanRegionsFile:
	input:
		coord="illumina/support_files/Multi-EthnicGlobal_D2.csv_Physical-and-Genetic-Coordinates.txt",
	output:
		regionsFile="data/snp_array/vcf/D2_cleanRegionsFile.txt",
	benchmark:
		repeat("benchmarks/cleanRegionsFile.array", 1)
	shell:
		"awk '{{print $2,$3}}' {input.coord} | "
		"sed 's/ /-/g' - | " # might be '\t' for GRCh37
		"sort -n - | "
		"sed '/0-0/d' - | sed '/XY/d' - | sed '/MT/d' - | sed '/Chr-MapInfo/d' | "
		"sed 's/-/\t/g' - > {output.regionsFile}"

rule extract_clean_regions:
	input:
		regionsFile="data/snp_array/vcf/D2_cleanRegionsFile.txt",
		vcf="data/snp_array/vcf/1-vcf-from-gtc.vcf.gz",
		vcfTBI="data/snp_array/vcf/1-vcf-from-gtc.vcf.gz.tbi",
		fasta=rules.bwa_index.output,
	output:
		vcf2="data/snp_array/vcf/2-vcf-from-cleanRegionsFile.vcf",
	benchmark:
		repeat("benchmarks/extract_clean_regions.array", 1)
	script:
		"../scripts/extract-cleanRegionsFile-for-2-vcf.py"

checkpoint download_lociNames_conversion_file:
	output:
		lociNamesConversionFile=protected("illumina/support_files/Multi-EthnicGlobal_D2_b150_rsids.txt"),
	params:
		rsidsUrl=config["urlSupportFiles"]["rsidConversion"],
		rzip=config["urlSupportFiles"]["rzip"],
		supportDir="illumina/support_files/",
	benchmark:
		repeat("benchmarks/download_lociNames_conversion_file.array", 1)
	shell:
		"""
		wget -t 0 --retry-connrefused -c -nH -P {params.supportDir} {params.rsidsUrl}
		sleep 30s
		unzip -o {params.supportDir}{params.rzip} -d {params.supportDir}
		sleep 30s
		rm {params.supportDir}{params.rzip}
		"""

rule replace_lociNames_to_rsIDs:
	input:
		vcf2="data/snp_array/vcf/2-vcf-from-cleanRegionsFile.vcf",
		lociNamesConversionFile="illumina/support_files/Multi-EthnicGlobal_D2_b150_rsids.txt",
	output:
		rsidsConversionFile="data/snp_array/vcf/B-rsids-from-Illumina-to-replace-locinames.txt",
	params:
		rsidsUrl=config["urlSupportFiles"]["rsidConversion"],
		rzip=config["urlSupportFiles"]["rzip"],
		rfile=config["urlSupportFiles"]["rfile"],
		workdir=config["workdirPath"],
	conda: "../../env/replace-lociNames-with-rsIDs.yaml",
	benchmark:
		repeat("benchmarks/replace_lociNames_to_rsIDs.array", 1)
	script:
		"../scripts/replace-lociNames-with-rsIDs.py"

rule original_lociNames:
	input:
		vcf2="data/snp_array/vcf/2-vcf-from-cleanRegionsFile.vcf",
	output:
		OriginalLociNames="data/snp_array/vcf/C-original-LociNames.txt",
	benchmark:
		repeat("benchmarks/original_lociNames.array", 1)
	shell:
		"awk '/^[^#]/ {{print $1,$2,$3}}' {input.vcf2} > {output.OriginalLociNames}"

rule insert_rsIDS:
	input:
		OriginalLociNames="data/snp_array/vcf/C-original-LociNames.txt",
		rsidsConversionFile="data/snp_array/vcf/B-rsids-from-Illumina-to-replace-locinames.txt",
		vcf2="data/snp_array/vcf/2-vcf-from-cleanRegionsFile.vcf",
	output:
		finalConversionFile="data/snp_array/vcf/D-chrom-pos-lociname-withPasted-matching-rsids.txt",
		vcf3="data/snp_array/vcf/3-lociNames-replaced-by-rsIDs.vcf",
		vcf3check="data/snp_array/vcf/3-vcf-check-for-correct-index.txt",
	benchmark:
		repeat("benchmarks/insert_rsIDS.array", 1)
	shell:
		"""
		paste {input.OriginalLociNames} {input.rsidsConversionFile} > {output.finalConversionFile}
		awk 'NR==FNR {{A[$1,$2] = $4; next}} ($1,$2) in A{{$3=A[$1,$2]}}1' OFS='\t' {output.finalConversionFile} {input.vcf2} > {output.vcf3}
		awk '/^[^#]/ {{print $1,$2,$3}}' {output.vcf3} > {output.vcf3check} # make sure rsids were correctly indexed
		"""

rule filter_for_dbSNP_vars:
	input:
		vcf3="data/snp_array/vcf/3-lociNames-replaced-by-rsIDs.vcf",
		rsidList="data/snp_array/vcf/E-rsids-found-from-dbSNP-all-file.txt",
	output:
		vcf4="data/snp_array/vcf/4-extracted-dbSNP-rsids.vcf",
	params:
		workdirPath=config['workdirPath'],
	benchmark:
		repeat("benchmarks/filter_for_dbSNP_vars.array", 1),
	resources:
		mem_mb=59900,
	script:
		"../scripts/filter-for-dbSNP-variants.py"