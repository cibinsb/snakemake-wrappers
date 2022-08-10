__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import uuid
import shlex
from pathlib import Path, PurePath
from snakemake.shell import shell

nf_version = snakemake.params.get("nf_version", "22.04.0-5697")
revision = snakemake.params.get("revision")
profile = snakemake.params.get("profile", [])
configs = snakemake.input.get("config", [])
with_tower = snakemake.params.get("with_tower")
extra = snakemake.params.get("extra", "")
work_dir = snakemake.params.get("w", None)
resume = snakemake.params.get("resume", False)

# placeholder value for custom work dir path
custom_work_dir = str(uuid.uuid4())
for output_file in str(snakemake.output).split():
    custom_work_dir = output_file.replace(".done", "")
    break
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
if work_dir:
    args += ["-w", work_dir]
if resume:
    args += ["-resume"]
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
single_dash_params = ["pipeline", "nf_version", "revision", "profile", "with_tower", "extra", "w", "resume"]
for name, value in snakemake.params.items():
    if name not in single_dash_params:
        add_parameter(name, value)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
args = " ".join(args)
pipeline = snakemake.params.pipeline

work_dir = str(Path(os.getcwd()) / Path("work"))
run_dir = PurePath(os.getcwd())
nfl_tmp_work = os.getenv("TMP_DIR")
nf_work = shlex.quote(str(Path(nfl_tmp_work) / run_dir.parent.name
                          / run_dir.name / Path("work") / Path(custom_work_dir)))
shell(
    """
    set -ue
    # setting NXF_WORK variable on the shell
    export NXF_WORK="{nf_work}"
    # Adding symbolic link to Nextflow work directory
    umask 0077
    mkdir -p "{nf_work}"
    mkdir -p "{work_dir}"
    ln -sf {nf_work} {work_dir}
    # Run nextflow
    #Same as "module load nextflow/{nf_version}", in case "module" is not available
    TOL_NEXTFLOW_DIR=/software/treeoflife/custom-installs/bin/nextflow/{nf_version}
    env NXF_JAVA_HOME=$TOL_NEXTFLOW_DIR/jdk-11 \
    $TOL_NEXTFLOW_DIR/nextflow/nextflow run {pipeline} {args} {extra} {log}
    """
)
