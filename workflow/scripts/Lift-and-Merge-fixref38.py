from os import path
from snakemake.shell import shell

# Extract arguments.
dbsnpDir = snakemake.params.get("dbsnpDir", "")
dbsnpFile = snakemake.params.get("dbsnpFile", "")
genomeFastaDir = snakemake.params.get("genomeFastaDir", "")
genomeFastaFile = snakemake.params.get("genomeFastaFile", "")
project = snakemake.params.get("project_from_wc", "")
refAssemblyVersion = snakemake.params.get("refAssemblyVersion_from_wc", "")
valid38 = snakemake.params.get("valid38_from_wc", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "mkdir -p data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/Temp/{valid38}/"
)

shell(
    "bcftools +fixref {snakemake.input.renamedChrs38}"
    " -O u --"
    " -d -f {genomeFastaDir}{genomeFastaFile}"
    " -i {dbsnpDir}{dbsnpFile} |"
    " bcftools sort"
    " -m 8G"
    " -T data/vcf_Lift-and-Merge/{project}/{refAssemblyVersion}/step6B-clean-vcfIDs38/Temp/{valid38}/"
    " -O v"
    " -o {snakemake.output.FixRef38}"
)