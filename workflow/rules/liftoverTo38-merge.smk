# 38
rule force_same_chr_scheme:
	input:
		vcfInput="data/vcf_merge/{vcf}.vcf",
	output:
		YesChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step1-chrConvention/{vcf}.ChrString.vcf",
	params:
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	conda: "../../env/bcftools1-14.yaml",
	benchmark:
		repeat("benchmarks/force_same_chr_scheme-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	script:
		"../scripts/force-same-chr-scheme-vcf.py"

rule extract_rsids:
	input:
		YesChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step1-chrConvention/{vcf}.ChrString.vcf",
	output:
		rsidList="data/vcf_merge/{project}/{refAssemblyVersion}/step2-rsids/{vcf}-rsidList.txt",
	params:
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	benchmark:
		repeat("benchmarks/extract_rsids-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		awk '/^[^#]/ {{print $3}}' {input.YesChrString} > {output.rsidList}
		"""

rule get_rsids_from_dbSNP:
	input:
		rsidList="data/vcf_merge/{project}/{refAssemblyVersion}/step2-rsids/{vcf}-rsidList.txt",
		tmpoFile="dbSNP/tempFile.to.remove",
	output:
		dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-rsids/{vcf}-dbSNP-IDs-from-rsIDList.txt",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbSNP']['file'],
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	shell:
		"""
		zgrep -w -f {input.rsidList} {params.dbsnpDir}{params.dbsnpFile} > {output.dbsnpExtractedIDsFile}
		"""

rule create_guide_file:
	input:
		dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-rsids/{vcf}-dbSNP-IDs-from-rsIDList.txt",
	output:
		guideFile="data/vcf_merge/{project}/{refAssemblyVersion}/step3-guideFile/{vcf}-guideFile-for-lift.txt",
	params:
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	script:
		"../scripts/guide-file-38-merge.py"

rule perform_liftover_to_38:
	input:
		guideFile="data/vcf_merge/{project}/{refAssemblyVersion}/step3-guideFile/{vcf}-guideFile-for-lift.txt",
		YesChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step1-chrConvention/{vcf}.ChrString.vcf",
	output:
		liftedOver="data/vcf_merge/{project}/{refAssemblyVersion}/step4-liftedOver/{vcf}.Lifted.vcf",
	params:
		refAssemblyVersion=config['Liftover']['desiredVersion'],
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	script:
		"../scripts/perform-lift-over-38.py"

rule zip:
	input:
		liftedOverVCF="data/vcf_merge/{project}/{refAssemblyVersion}/step4-liftedOver/{vcf}.Lifted.vcf",
	output:
		liftedOverGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step4-liftedOver/{vcf}.Lifted.vcf.gz",
	params:
	resources:
		mem_mb=1500,
		runtime="01:00:00",
	shell:
		"""
		bcftools view -O z -o {output.liftedOverGZ} {input.liftedOverVCF}
		"""

rule index:
	input:
		liftedOverGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step4-liftedOver/{vcf}.Lifted.vcf.gz",
	output:
		liftedOverGZ_TBI="data/vcf_merge/{project}/{refAssemblyVersion}/step4-liftedOver/{vcf}.Lifted.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="01:00:00",
	shell:
		"""
		bcftools index --tbi {input.liftedOverGZ}
		"""
