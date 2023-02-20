from os import path
from snakemake.shell import shell

# Extract arguments.
dbsnpDir = snakemake.params.get("dbsnpDir", "")
dbsnpFile = snakemake.params.get("dbsnpFile", "")
genomeFastaDir = snakemake.params.get("genomeFastaDir", "")
genomeFastaFile = snakemake.params.get("genomeFastaFile", "")
project = snakemake.params.get("project_from_wc", "")
refAssemblyVersion = snakemake.params.get("refAssemblyVersion_from_wc", "")
vcf37 = snakemake.params.get("vcf37_from_wc", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "mkdir -p data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/Temp/{vcf37}/"
)

shell(
    "bcftools +fixref {snakemake.input.cleanDBSNPmyData}"
    " -O u --"
    " -d -f {genomeFastaDir}{genomeFastaFile}"
    " -i {dbsnpDir}{dbsnpFile} |"
    " bcftools sort"
    " -m 8G"
    " -T data/vcf_Merge-and-Lift/{project}/{refAssemblyVersion}/step6A-clean-vcfIDs37/Temp/{vcf37}/"
    " -O v"
    " -o {snakemake.output.FixRef}"
)