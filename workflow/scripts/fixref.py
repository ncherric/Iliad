from os import path
from snakemake.shell import shell

# Extract arguments.
dbsnpDir=snakemake.params.get("dbsnpDir", "")
dbsnpFile=snakemake.params.get("dbsnpFile", "")
project=snakemake.params.get("project_from_wc", "")
refAssemblyVersion=snakemake.params.get("refAssemblyVersion_from_wc", "")
vcf=snakemake.params.get("vcf_from_wc", "")

reference = snakemake.input.ref

# if isinstance(reference, str):
# 	reference = path.splitext(snakemake.input.ref)[0]
# else:
# 	reference = path.splitext(snakemake.input.ref[0])[0]

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "mkdir -p data/vcf_merge/{project}/{refAssemblyVersion}/step7-fixref/Temp/{vcf}/"
)

shell(
    "bcftools +fixref {snakemake.input.cleanDBSNPmyData}"
    " -O u --"
    " -d -f {reference}"
    # " -i {snakemake.input.dbsnpExtractedIDsFile} |"
    " -i {dbsnpDir}{dbsnpFile} |"
    " bcftools sort"
    " -m 8G"
    " -T data/vcf_merge/{project}/{refAssemblyVersion}/step7-fixref/Temp/{vcf}/"
    " -O v"
    " -o {snakemake.output.FixRef}"
)