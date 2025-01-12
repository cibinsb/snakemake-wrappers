__author__ = "Christopher Schröder, Patrik Smeds"
__copyright__ = "Copyright 2022, Christopher Schröder, Patrik Smeds"
__email__ = "christopher.schroeder@tu-dortmund.de, patrik.smeds@gmail.com"
__license__ = "MIT"

from os import path

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Check inputs/arguments.
if len(snakemake.input) == 0:
    raise ValueError("A reference genome has to be provided.")
elif len(snakemake.input) > 1:
    raise ValueError("Please provide exactly one reference genome as input.")

valid_suffixes = {
    ".0123",
    ".amb",
    ".ann",
    ".pac",
    ".pos_packed",
    ".suffixarray_uint64",
    ".suffixarray_uint64_L0_PARAMETERS",
    ".suffixarray_uint64_L1_PARAMETERS",
    ".suffixarray_uint64_L2_PARAMETERS",
}


def get_valid_suffix(path):
    for suffix in valid_suffixes:
        if path.endswith(suffix):
            return suffix


prefixes = set()
for s in snakemake.output:
    suffix = get_valid_suffix(s)
    if suffix is None:
        raise ValueError(f"{s} cannot be generated by bwa-meme index (invalid suffix).")
    prefixes.add(s[: -len(suffix)])

if len(prefixes) != 1:
    raise ValueError("Output files must share common prefix up to their file endings.")
(prefix,) = prefixes

shell(
    "(bwa-meme index -a meme -p {prefix} {snakemake.input[0]} -t {snakemake.threads} && build_rmis_dna.sh {snakemake.input[0]}) {log}"
)
