from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

out_args = ""

out_file = snakemake.output.get("out")
if out_file:
    out_args += f" -o {out_file}"
elif len(snakemake.output) == 1:
    out_args += f" -o {snakemake.output[0]}"

faa_file = snakemake.output.get("faa")
if faa_file:
    out_args += f" -a {faa_file}"

fna_file = snakemake.output.get("fna")
if fna_file:
    out_args += f" -d {fna_file}"

stat_file = snakemake.output.get("stat")
if stat_file:
    out_args += f" -s {stat_file}"

shell("prodigal " "-i {snakemake.input.fasta} " "{out_args} " "{extra} " "{log}")
