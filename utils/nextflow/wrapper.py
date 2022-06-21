__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import shlex
from pathlib import Path, PurePath
from snakemake.shell import shell

revision = snakemake.params.get("revision")
profile = snakemake.params.get("profile", [])
configs = snakemake.params.get("config", [])
with_tower = snakemake.params.get("with_tower")
extra = snakemake.params.get("extra", "")
if isinstance(profile, str):
    profile = [profile]

args = []

if revision:
    args += ["-revision", revision]
if profile:
    args += ["-profile", ",".join(profile)]
if configs:
    for config in configs:
        args.extend(["-config", config])
if with_tower:
    args += ["-with-tower", with_tower]
print(args)

# TODO pass threads in case of single job
# TODO limit parallelism in case of pipeline
# TODO handle other resources

add_parameter = lambda name, value: args.append("--{} {}".format(name, value))

for name, files in snakemake.input.items():
    if isinstance(files, list):
        # TODO check how multiple input files under a single arg are usually passed to nextflow
        files = ",".join(files)
    add_parameter(name, files)
single_dash_params = ["pipeline", "revision", "profile", "config", "with_tower", "extra"]
for name, value in snakemake.params.items():
    if name not in single_dash_params:
        add_parameter(name, value)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
args = " ".join(args)
pipeline = snakemake.params.pipeline

work_dir = os.getcwd()
run_dir = PurePath(work_dir)
nfl_tmp_work = os.getenv("TMP_DIR")
nf_work = shlex.quote(str(Path(nfl_tmp_work) / run_dir.parent.name / run_dir.name / Path("work")))

# setting NXF_WORK variable on the shell
run_command = f'export NXF_WORK="{nf_work}"; '
# Adding symbolic link to Nextflow work directory
run_command += f'set -ue; umask 0077; mkdir -p {nf_work}; \
                mkdir -p {work_dir}; \
                ln -sf {nf_work} {work_dir}'
shell(
    """
    {run_command} | nextflow run {pipeline} {args} {extra} {log}
    """
)
