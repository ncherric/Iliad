from snakemake.shell import shell

# Extract arguments.
# extra = snakemake.params.get("extra", "")
refAssemblyVersion = snakemake.params.get("refAssemblyVersion", ""),

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
	"awk 'NR==FNR {{A[$1,$3] = $2; next}} ($1,$3) in A{{$2=A[$1,$3]}}1' OFS='\t' {snakemake.input.guideFile} {snakemake.input.YesChrString} > {snakemake.output.liftedOver}"
)
