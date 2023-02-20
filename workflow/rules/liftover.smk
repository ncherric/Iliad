# ruleorder: unzip > query1_mergeData_ID_list > combine_MyData_SNPs > get_rsids_from_dbSNP > rsids_from_dbSNP_index > query2_mergeData_ID_list > get_clean_rsids_from_MyData > clean_rsids_for_MyData_index > force_same_chr_scheme
ruleorder: zip1 > get_rsids_from_dbSNP > perform_liftover > sort > fix_ref > merge_vcfs

rule zip1: # change this to zip
	input:
		unzipped="data/vcf_merge/{vcf}.vcf",
	output:
		zipped="data/vcf_merge/{vcf}.zipped.vcf.gz",
	params:
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	# conda: "../../env/bcftools1-14.yaml",
	benchmark:
		repeat("benchmarks/zip1-{vcf}.merger", 1)
	shell:
		"""
		bcftools view -O z -o {output.zipped} {input.unzipped}
		"""

rule index1:
	input:
		zipped1="data/vcf_merge/{vcf}.zipped.vcf.gz",
	output:
		indexed1="data/vcf_merge/{vcf}.zipped.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/index1-{vcf}.merger", 1)
	shell:
		"""
		bcftools index --tbi {input.zipped1}
		"""

rule query1_mergeData_ID_list:
	input:
		zipped1="data/vcf_merge/{vcf}.zipped.vcf.gz",
		indexed1="data/vcf_merge/{vcf}.zipped.vcf.gz.tbi",
	output:
		rsidList="data/vcf_merge/{project}/{refAssemblyVersion}/step1-MyData-SNP-IDs/one.{vcf}.rsIDs.txt",
		# vcfGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step1-MyData-SNP-IDs/{vcf}.vcf.gz",
		# vcfTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step1-MyData-SNP-IDs/{vcf}.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query1_mergeData_ID_list-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	script:
		"../scripts/query1-mergeData-ID-list.py"
	# shell:
	# 	"""
	# 	bcftools view ^-t chr0 -O u {input.vcfInput} | \\
	# 	bcftools query -f %ID\\\\n > {output.rsidList}
	# 	bcftools view -O z -o {output.vcfGZ} {input.vcfInput}
	# 	bcftools index --tbi {output.vcfGZ}
	# 	"""

rule combine_MyData_SNPs:
	input:
		rsidList=get_combineMyData_filepaths,
	output:
		combinedSNPlist="data/vcf_merge/{project}/{refAssemblyVersion}/step1-MyData-SNP-IDs/combinedMyData.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/combine_MyData_SNPs-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		cat {input} | sort -V - | uniq - > {output.combinedSNPlist}
		"""
	# script:
	# 	"../scripts/combine-MyData-SNPs.py"

rule get_rsids_from_dbSNP:
	input:
		combinedSNPlist="data/vcf_merge/{project}/{refAssemblyVersion}/step1-MyData-SNP-IDs/combinedMyData.rsIDs.txt",
	output:
		dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbSNP']['file'], # change to 37
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/get_rsids_from_dbSNP-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools view -i ID=@{input.combinedSNPlist} -O z -o {output.dbsnpExtractedIDsFile} {params.dbsnpDir}{params.dbsnpFile}
		"""

rule rsids_from_dbSNP_index:
	input:
		dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz",
	output:
		dbsnpExtractedIDsFileTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/rsids_from_dbSNP_index-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools index --tbi {input.dbsnpExtractedIDsFile}
		"""

rule query2_mergeData_ID_list:
	input:
		dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz",
		dbsnpExtractedIDsFileTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz.tbi",
	output:
		cleanDBSNPlist="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/two.clean-dbSNP-combinedMyData.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query2_mergeData_ID_list-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.dbsnpExtractedIDsFile} > {output.cleanDBSNPlist}
		"""

rule get_clean_rsids_from_MyData:
	input:
		cleanDBSNPlist="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/two.clean-dbSNP-combinedMyData.rsIDs.txt",
		vcfGZ="data/vcf_merge/{vcf}.zipped.vcf.gz",
		vcfTBI="data/vcf_merge/{vcf}.zipped.vcf.gz.tbi",
	output:
		cleanDBSNPmyData="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/get_clean_rsids_from_MyData-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools view -i ID=@{input.cleanDBSNPlist} -O z -o {output.cleanDBSNPmyData} {input.vcfGZ}
		"""

rule clean_rsids_for_MyData_index:
	input:
		cleanDBSNPmyData="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz",
	output:
		cleanDBSNPmyDataTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="12:00:00",
	benchmark:
		repeat("benchmarks/clean_rsids_for_MyData_index-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools index --tbi {input.cleanDBSNPmyData}
		"""

########################################## just the VCF for adding chr
# rule get_clean_rsids_from_MyData:
# 	input:
# 		cleanDBSNPlist="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/two.clean-dbSNP-combinedMyData.rsIDs.txt",
# 		vcfGZ="data/vcf_merge/{vcf}.zipped.vcf.gz",
# 		vcfTBI="data/vcf_merge/{vcf}.zipped.vcf.gz.tbi",
# 	output:
# 		cleanDBSNPmyData="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz",
# 	resources:
# 		mem_mb=1500,
# 		runtime="12:00:00",
# 	benchmark:
# 		repeat("benchmarks/get_clean_rsids_from_MyData-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools view -i ID=@{input.cleanDBSNPlist} -O z -o {output.cleanDBSNPmyData} {input.vcfGZ}
# 		"""






# rule force_same_chr_scheme:
# 	input:
# 		cleanDBSNPmyData="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz",
# 		cleanDBSNPmyDataTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz.tbi",
# 		renameChrFile="config/renameChrFile-2.txt",
# 	output:
# 		NoChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step4-chrConvention/ChrStringNone.{vcf}.vcf",
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:30:00",
# 	benchmark:
# 		repeat("benchmarks/force_same_chr_scheme-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools annotate --rename-chrs {input.renameChrFile} {input.cleanDBSNPmyData} -O v -o {output.NoChrString}
# 		"""

###### #######3 ####### ###### keep this one
# rule force_same_chr_scheme:
# 	input:
# 		cleanDBSNPmyData="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz",
# 		cleanDBSNPmyDataTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz.tbi",
# 	output:
# 		ChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step4-chrConvention/{vcf}.ChrString.vcf",
# 	params:
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:30:00",
# 	benchmark:
# 		repeat("benchmarks/force_same_chr_schemeNoah-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	script:
# 		"../scripts/force-same-chr-scheme-vcf.py"


# rule force_same_chr_scheme:
# 	input:
# 		vcfInput="data/vcf_merge/{vcf}.vcf",
# 	output:
# 		NoChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step1-chrConvention/ChrStringNone.{vcf}.vcf",
# 	params:
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:30:00",
# 	conda: "../../env/bcftools1-14.yaml",
# 	benchmark:
# 		repeat("benchmarks/force_same_chr_scheme-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	script:
# 		"../scripts/force-same-chrNone-scheme-vcf.py"

# rule extract_rsids:
# 	input:
# 		NoChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step1-chrConvention/ChrStringNone.{vcf}.vcf",
# 	output:
# 		rsidList="data/vcf_merge/{project}/{refAssemblyVersion}/step2-rsids/{vcf}-rsidList.txt",
# 	params:
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:30:00",
# 	benchmark:
# 		repeat("benchmarks/extract_rsids-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		awk '/^[^#]/ {{print $3}}' {input.NoChrString} > {output.rsidList}
# 		"""

# rule get_rsids_from_dbSNP:
# 	input:
# 		rsidList="data/vcf_merge/{project}/{refAssemblyVersion}/step2-rsids/{vcf}-rsidList.txt",
# 		tmpoFile="dbSNP/tempFile.to.remove",
# 	output:
# 		dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-rsids/{vcf}-dbSNP-IDs-from-rsIDList.txt",
# 	params:
# 		dbsnpDir="dbSNP/",
# 		dbsnpFile=config['dbSNP']['file'],
# 	resources:
# 		mem_mb=1500,
# 		runtime="06:00:00",
# 	shell:
# 		"""
# 		zgrep -w -f {input.rsidList} {params.dbsnpDir}{params.dbsnpFile} > {output.dbsnpExtractedIDsFile}
# 		"""



# KEEPING ABOVE FOR NOW

###################################################################

# NEED BELOW


# rule zip2:
# 	input:
# 		ChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step4-chrConvention/{vcf}.ChrString.vcf",
# 	output:
# 		ChrStringz="data/vcf_merge/{project}/{refAssemblyVersion}/step4-chrConvention/{vcf}.ChrString.vcf.gz",
# 	resources:
# 		mem_mb=1500,
# 		runtime="06:00:00",
# 	benchmark:
# 		repeat("benchmarks/zip2-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools view -O z -o {output.ChrStringz} {input.ChrString}
# 		"""

# rule index2:
# 	input:
# 		ChrStringz="data/vcf_merge/{project}/{refAssemblyVersion}/step4-chrConvention/{vcf}.ChrString.vcf.gz",
# 	output:
# 		ChrStringt="data/vcf_merge/{project}/{refAssemblyVersion}/step4-chrConvention/{vcf}.ChrString.vcf.gz.tbi",
# 	resources:
# 		mem_mb=1500,
# 		runtime="06:00:00",
# 	benchmark:
# 		repeat("benchmarks/index2-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools index --tbi {input.ChrStringz}
# 		"""

rule fix_ref:
	input:
		cleanDBSNPmyData="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz",
		cleanDBSNPmyDataTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.vcf.gz.tbi",
		# sortedGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
		# sortedTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz.csi",
		# ref=rules.bwa_index.output,
		ref="resources/human_g1k_v37.fasta",
		tmpoFile="dbSNP/tempFile.to.remove",
		# dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz",
		# dbsnpExtractedIDsFileTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz.tbi",
	output:
		FixRef="data/vcf_merge/{project}/{refAssemblyVersion}/step4-fixref/fixref.{vcf}.vcf",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbSNP']['file'],
		project_from_wc=lambda wc: wc.get("project"),
		refAssemblyVersion_from_wc=lambda wc: wc.get("refAssemblyVersion"),
		vcf_from_wc=lambda wc: wc.get("vcf"),
	resources:
		mem_mb=9000,
		runtime="01:00:00",
	benchmark:
		repeat("benchmarks/fix_ref-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	script:
		"../scripts/fixref.py"

# rule fix_ref_index:
# 	input:
# 		FixRefGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step4-fixref/fixref.{vcf}.vcf.gz",
# 	output:
# 		FixRefTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step4-fixref/fixref.{vcf}.vcf.gz.tbi",
# 	resources:
# 		mem_mb=1500,
# 		runtime="00:30:00",
# 	benchmark:
# 		repeat("benchmarks/fix_ref_index-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools index --tbi {input.FixRefGZ} -o {output.FixRefTBI}
# 		"""



rule create_guide_file:
	input:
		dbsnpExtractedIDsFile="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz",
		dbsnpExtractedIDsFileTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step2-dbSNP-IDs/one.dbSNP-combinedMyData.vcf.gz.tbi",
		# cleanDBSNPmyDataZ="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.zipped.vcf.gz",
		# cleanDBSNPmyDataT="data/vcf_merge/{project}/{refAssemblyVersion}/step3-dbSNP-IDs/one.{vcf}.cleanSNPs.zipped.vcf.gz.tbi",
	output:
		guideFile="data/vcf_merge/{project}/{refAssemblyVersion}/step5-guideFile/guideFile-for-lift.txt",
	params:
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/create_guide_file-{project}-{refAssemblyVersion}.merger", 1)
	script:
		"../scripts/guide-file-merge.py"

rule perform_liftover:
	input:
		guideFile="data/vcf_merge/{project}/{refAssemblyVersion}/step5-guideFile/guideFile-for-lift.txt",
		# ChrString="data/vcf_merge/{project}/{refAssemblyVersion}/step4-chrConvention/{vcf}.ChrString.vcf",
		FixRef="data/vcf_merge/{project}/{refAssemblyVersion}/step4-fixref/fixref.{vcf}.vcf",
	output:
		liftedOver="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/one.Lifted.{vcf}.vcf",
	params:
		refAssemblyVersion=config['Liftover']['desiredVersion'],
	resources:
		mem_mb=1500,
		runtime="06:00:00",
	benchmark:
		repeat("benchmarks/perform_liftover-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	script:
		"../scripts/perform-lift-over.py"

rule sort:
	input:
		liftedOverVCF="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/one.Lifted.{vcf}.vcf",
	output:
		sortedGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
	params:
		tempDir="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/TempDir/{vcf}/"
	resources:
		mem_mb=4000,
		runtime="01:00:00",
	benchmark:
		repeat("benchmarks/zip-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		mkdir -p {params.tempDir}
		bcftools sort -m 3G -T {params.tempDir} -O z -o {output.sortedGZ} {input.liftedOverVCF}
		"""

rule sort_index:
	input:
		sortedGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
	output:
		sortedTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz.csi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:40:00",
	benchmark:
		repeat("benchmarks/index-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools index --csi {input.sortedGZ}
		"""


rule filter:
	input:
		sortedGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
		sortedTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz.csi",
	output:
		filtered="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/three.Filtered.{vcf}.vcf.gz",
	params:
	resources:
		mem_mb=1500,
		runtime="00:40:00",
	benchmark:
		repeat("benchmarks/filter-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools view -i 'F_MISSING < 0.05' {input.sortedGZ} -O z -o {output.filtered}
		"""

rule filter_index:
	input:
		filtered="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/three.Filtered.{vcf}.vcf.gz",
	output:
		filteredTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/three.Filtered.{vcf}.vcf.gz.tbi",
	params:
	resources:
		mem_mb=1500,
		runtime="00:40:00",
	benchmark:
		repeat("benchmarks/filter_index-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
	shell:
		"""
		bcftools index --tbi {input.filtered}
		"""


# rule zip:
# 	input:
# 		sorted="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
# 	output:
# 		sortedGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
# 	params:
# 	resources:
# 		mem_mb=1500,
# 		runtime="01:00:00",
# 	benchmark:
# 		repeat("benchmarks/zip-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools view -O z -o {output.sortedGZ} {input.sorted}
# 		"""

# rule index:
# 	input:
# 		sortedGZ="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz",
# 	output:
# 		sortedTBI="data/vcf_merge/{project}/{refAssemblyVersion}/step6-liftedOver/two.Sorted.{vcf}.vcf.gz.tbi",
# 	params:
# 	resources:
# 		mem_mb=1500,
# 		runtime="01:00:00",
# 	benchmark:
# 		repeat("benchmarks/index-{vcf}-{project}-{refAssemblyVersion}.merger", 1)
# 	shell:
# 		"""
# 		bcftools index --tbi {input.sortedGZ}
# 		"""
