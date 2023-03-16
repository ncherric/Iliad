from glob import glob
from pathlib import Path

def aggregate_38_VCFs(wildcards):
	print(f"aggregate_38_VCFs is: {wildcards}")
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	print(f"this is checkpoint_output: {checkpoint_output}")
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	print(f"ids is: {ids}")
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	print(f"baseNames is: {baseNames}")
	files = expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step3A-InputVCFs-38/{{vcf38}}.vcf.gz", vcf38=baseNames)
	print(f"files is: {files}")
	return files

rule annotate_38_VCFs:
	input:
		file38=aggregate_38_VCFs,
		dbsnpTemp="dbSNP/tempFile38.to.remove",
	output:
		annotated38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/annotated.{vcf38}.vcf.gz",
		annotated38index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/annotated.{vcf38}.vcf.gz.tbi",
	params:
		dbsnpDir="dbSNP/",
		dbsnpFile=config['dbsnpLiftMerge']['file38'],
	resources:
		mem_mb=3500,
		runtime="04:30:00",
	benchmark:
		repeat("benchmarks/annotate_38_VCFs-{vcf38}-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		bcftools annotate -a {params.dbsnpDir}{params.dbsnpFile} -c ID -O z -o {output.annotated38} {input.file38}
		bcftools index --tbi {output.annotated38}
		"""

rule query_38_VCF_IDs:
	input:
		annotated38="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/annotated.{vcf38}.vcf.gz",
		annotated38index="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/annotated.{vcf38}.vcf.gz.tbi",
	output:
		rsidList="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/queried.{vcf38}.rsIDs.txt"
	resources:
		mem_mb=1500,
		runtime="00:05:00",
	benchmark:
		repeat("benchmarks/query_38_VCF_IDs-{vcf38}-{{project}}-{{refAssemblyVersion}}.LiftAndMerge", 1)
	shell:
		"""
		bcftools query -f %ID\\\\n {input.annotated38} > {output.rsidList}
		"""

def get_combineMyData_38filepaths(wildcards):
	checkpoint_output = os.path.dirname(checkpoints.check_basenames_v38.get(**wildcards).output[0])
	ids = glob(join(checkpoint_output, "38path_*.vcf.gz"))
	baseNames = []
	for id in ids:
		baseNames.append(Path(id).stem.rsplit('.',maxsplit=1)[0])
	return expand("data/vcf_Lift-and-Merge/{{project}}/{{refAssemblyVersion}}/step4A-Input-vcfIDs38/queried.{vcf38}.rsIDs.txt",vcf38=baseNames)


rule combine_MyData_SNPs:
	input:
		rsidList=get_combineMyData_38filepaths,
	output:
		combinedSNPlist="data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step4A-Input-vcfIDs38/combinedMyData.rsIDs.txt",
	resources:
		mem_mb=1500,
		runtime="00:10:00",
	benchmark:
		repeat("benchmarks/combine_MyData_SNPs-{project}-{refAssemblyVersion}.LiftAndMerge", 1)
	shell:
		"""
		cat {input} | sort -V - | uniq - > {output.combinedSNPlist}
		"""