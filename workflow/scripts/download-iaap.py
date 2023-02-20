from os import path
import re
import tempfile
from snakemake.shell import shell


# Extract arguments.
# extra = snakemake.params.get("extra", "")
productDir=snakemake.params.get("productDir", "")
bpm=snakemake.params.get("bpm", "")
mzip=snakemake.params.get("mzip", "")
egt=snakemake.params.get("egt", "")
czip=snakemake.params.get("czip", "")

iaapIlluminaDownloadPath=snakemake.params.get("iaapIlluminaDownloadPath", "")
iaapIlluminaDownloadTar=snakemake.params.get("iaapIlluminaDownloadTar", "")
iaapIlluminaDownload=snakemake.params.get("iaapIlluminaDownload", "")
iaapcli=snakemake.params.get("iaapcli", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    	"mkdir -p {snakemake.output.iaapDir}"
)

shell(
		"wget -t 0 --retry-connrefused -c -nH {iaapIlluminaDownloadPath} -P {snakemake.output.iaapDir}"
)

shell(
		"tar -xzvf illumina/gencall/{iaapIlluminaDownloadTar} -C illumina/gencall/ {iaapIlluminaDownload}/{iaapcli} --strip-components=1"
)

shell(
		"rm illumina/gencall/{iaapIlluminaDownloadTar}"
)

shell(
		"chmod 777 {snakemake.output.iaap}"
)
