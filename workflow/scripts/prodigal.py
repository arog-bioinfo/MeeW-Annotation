from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

out_file = snakemake.output.get("out") or ""
if not out_file and len(snakemake.output) == 1:
    out_file = snakemake.output[0]
out_option = "-o" if out_file else ""

faa_file = snakemake.output.get("faa") or ""
faa_option = "-a" if faa_file else ""

fna_file = snakemake.output.get("fna") or ""
fna_option = "-d" if fna_file else ""

stat_file = snakemake.output.get("stat") or ""
stat_option = "-s" if stat_file else ""

shell(
    "prodigal "
    "-i {snakemake.input.fasta:q} "
    "{out_option} {out_file:q} "
    "{faa_option} {faa_file:q} "
    "{fna_option} {fna_file:q} "
    "{stat_option} {stat_file:q} "
    "{extra} "
    "{log}"
)
