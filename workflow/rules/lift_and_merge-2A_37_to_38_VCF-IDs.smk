from glob import glob
from pathlib import Path

def aggregate_37_VCFs(wildcards):
	print(f"aggregate_37_VCFs is: {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v37.get(**wildcards).output[0])
	print(f"this is checkpoint_output: {checkpoint_output}")
	ids = glob(join(checkpoint_output, "37path_*.vcf.gz"))
	print(f"ids is: {ids}")
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	print(f"baseNames is: {baseNames}")
	files = expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step3A-InputVCFs-37/{{vcf37}}.vcf.gz", vcf37=baseNames)
	print(f"files is: {files}")
	return files

rule annotate_37_VCFs:
	input:
		file37=aggregate_37_VCFs,
		dbsnpTemp="dbSNP/tempFile37.to.remove",
	output:
		annotated37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/annotated.{vcf37}.vcf.gz",
		annotated37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/annotated.{vcf37}.vcf.gz.tbi",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file37'],
	resources:
		mem_mb=3500,
		runtime="04:30:00",
	benchmark:
		repeat("benchmarks/annotate_37_VCFs-{vcf37}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate -a {params.dbsnpDir}{params.dbsnpFile} -c ID -O z -o {output.annotated37} {input.file37}
		bcftools index --tbi {output.annotated37}
		"""

rule query_37_VCF_IDs:
	input:
		annotated37="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/annotated.{vcf37}.vcf.gz",
		annotated37index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/annotated.{vcf37}.vcf.gz.tbi",
	output:
		rsidList="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/queried.{vcf37}.rsIDs.txt"
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query_37_VCF_IDs-{vcf37}-{{project}}-{{refAssemblyVersion}}.LiftAndMerge", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.annotated37} > {output.rsidList}
		"""

def get_combineMyData_37filepaths(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v37.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "37path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step4A-Input-vcfIDs37/queried.{vcf37}.rsIDs.txt",vcf37=baseNames)


rule combine_MyData_SNPs:
	input:
		rsidList=get_combineMyData_37filepaths,
	output:
		combinedSNPlist="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs37/combinedMyData.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/combine_MyData_SNPs-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input} | sort -V - | uniq - > {output.combinedSNPlist}
		"""