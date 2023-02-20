# Step 0 - retrieve necessary Illumina support files for your SNP array
checkpoint download_product_files:
	output:
		manifestFile="illumina/product_files/Multi-EthnicGlobal_D2.bpm",
		clusterFile="illumina/product_files/Multi-EthnicGlobal_D1_ClusterFile.egt",
	params:
		productDir="illumina/product_files/",
		bpm=config["urlProductFiles"]["manifest"],
		mzip=config["urlProductFiles"]["mzip"],
		egt=config["urlProductFiles"]["cluster"],
		czip=config["urlProductFiles"]["czip"],
	# benchmark:
	# 	repeat("benchmarks/download_product_files.array", 3)
	script:
		"../scripts/download-product-files.py"

# Step 1 - convert .idat files to .gtc files
# iaap-cli

checkpoint download_iaap:
	output:
		iaapDir=directory("illumina/gencall/"),
		iaap="illumina/gencall/iaap-cli/iaap-cli",
	params:
		iaapIlluminaDownloadPath=config['Illumina']['ftpDownload'],
		iaapIlluminaDownloadTar=config['Illumina']['DownloadTarFile'],
		iaapIlluminaDownload=config['Illumina']['Download'],
		iaapcli=config['Illumina']['iaapcli'],
	# benchmark:
	# 	repeat("benchmarks/download_iaap.array", 3)
	script:
		"../scripts/download-iaap.py"

rule gencall:
	input:
		bpm="illumina/product_files/Multi-EthnicGlobal_D2.bpm", # be sure that the latest manifest is used. manifest D1 was most recent in this build - D2 for GRCh38
		egt="illumina/product_files/Multi-EthnicGlobal_D1_ClusterFile.egt",
		idat="data/snp_array/idat/",
		iaapcli="illumina/gencall/iaap-cli/iaap-cli"
	output:
		gtc=directory("data/snp_array/gtc/"),
	# benchmark:
	# 	repeat("benchmarks/gencall.array", 3)
	shell:
		"""
		mkdir -p {output.gtc}
		{input.iaapcli} gencall {input.bpm} {input.egt} {output.gtc} -f {input.idat} -g
		"""

