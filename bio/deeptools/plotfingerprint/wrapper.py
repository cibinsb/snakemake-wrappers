__author__ = "Antonie Vietor"
__copyright__ = "Copyright 2020, Antonie Vietor"
__email__ = "antonie.v@gmx.de"
__license__ = "MIT"

import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

out_counts = snakemake.output.get("counts")

optional_output = ""

if out_counts:
    optional_output += " --outRawCounts {out_counts} ".format(out_counts=out_counts)

os.system(
    f"(plotFingerprint "
    "-b {snakemake.input.bam_files} "
    "-o {snakemake.output.fingerprint} "
    "{optional_output} "
    "{snakemake.params}) {log}"
)