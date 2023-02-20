from os.path import splitext

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Check inputs/arguments.
if len(snakemake.input) == 0:
    raise ValueError("A reference genome has to be provided!")
elif len(snakemake.input) > 1:
    raise ValueError("Only one reference genome can be inputed!")

# Prefix that should be used for the database
prefix = snakemake.params.get("prefix", splitext(snakemake.output[0])[0])

if len(prefix) > 0:
    prefix = "-p " + prefix

# Contrunction algorithm that will be used to build the database, default is bwtsw
construction_algorithm = snakemake.params.get("algorithm", "")

if len(construction_algorithm) != 0:
    construction_algorithm = "-a " + construction_algorithm

shell(
    "bwa index" " {prefix}" " {construction_algorithm}" " {snakemake.input[0]}" " {log}"
)