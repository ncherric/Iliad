from os import path
import re
import tempfile
from snakemake.shell import shell


# Extract arguments.
extra = snakemake.params.get("extra", "")

sort = snakemake.params.get("sorting", "none")
sort_order = snakemake.params.get("sort_order", "coordinate")
sort_extra = snakemake.params.get("sort_extra", "")

reference = snakemake.input.ref
if isinstance(reference, str):
	reference = path.splitext(snakemake.input.ref)[0]
else:
	reference = path.splitext(snakemake.input.ref[0])[0]


if re.search(r"-T\b", sort_extra) or re.search(r"--TMP_DIR\b", sort_extra):
	sys.exit(
		"You have specified temp dir (`-T` or `--TMP_DIR`) in params.sort_extra; this is automatically set from params.tmp_dir."
	)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


# Check inputs/arguments.
if not isinstance(snakemake.input.reads, str) and len(snakemake.input.reads) not in {
	1,
	2,
}:
	raise ValueError("input must have 1 (single-end) or " "2 (paired-end) elements")

pipe_cmd = (
	"/usr/bin/samtools-1.15/bin/samtools sort {sort_extra}"
	" -o {snakemake.output.sortedBam} -@ {snakemake.threads}"
)

with tempfile.TemporaryDirectory() as tmp:
	shell(
		"(bwa mem"
		" -t {snakemake.threads}"
		" {extra}"
		" {reference}"
		" {snakemake.input.reads}"
		" | " + pipe_cmd + ") {log}"
	)


shell(
	"/usr/bin/samtools-1.15/bin/samtools index -b -@ {snakemake.threads} {snakemake.output.sortedBam}"
) #########
#############
#############    FIX THIS 
#############
#############    INDEX WITH SAMTOOLS
