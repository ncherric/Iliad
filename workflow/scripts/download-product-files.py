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

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


shell(
		"mkdir -p {productDir}"
)

shell(
		"wget -t 0 --retry-connrefused -c -nH -P {productDir} {bpm}"
)

shell(
		"sleep 30s"
)

shell(
		"unzip -o {productDir}{mzip} -d {productDir}"
)

shell(
		"sleep 30s"
)

shell(
		"rm {productDir}{mzip}"
)

shell(
		"sleep 5s"
)

shell(
		"wget -t 0 --retry-connrefused -c -nH -P {productDir} {egt}"
)

shell(
		"sleep 30s"
)

shell(
		"unzip -o {productDir}{czip} -d {productDir}"
)

shell(
		"sleep 30s"
)

shell(
		"rm {productDir}{czip}"
)

shell(
		"sleep 5s"
)