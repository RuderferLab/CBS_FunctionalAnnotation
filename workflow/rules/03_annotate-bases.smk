from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["MAIN"]["IO_DIRS"]["PROCESS_DIRS"]["STEP_3"]

# Other
CBS_ACTIVITY = Path("results/02_process-motifs/MA0139.2-track.masked.pwm-activity.bed")

# Conservation files
LINSIGHT_SCORES = config["MAIN"]["CONSERVATION"]["LINSIGHT"]
PHSTCONS_SCORES = config["MAIN"]["CONSERVATION"]["PHASTCONS"]
GERP_SCORES = config["MAIN"]["CONSERVATION"]["GERP"]
PHYLOP_SCORES = config["MAIN"]["CONSERVATION"]["PHYLOP"]


# ------------- #
# Rules         #
# ------------- #

rule all:
    input:
        Path(PROCESS_DIR, "conservation.bed"),


# ------------- #
#     ----      #
# ------------- #


rule explode_input:
    """
    Explodes control track into positions
    """
    input:
        CBS_ACTIVITY,
    output:
        Path(
            PROCESS_DIR,
            "positions.bed",
        ),
    params:
        motif_length=15,
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        stdout="workflow/logs/explode_input.stdout",
        stderr="workflow/logs/explode_input.stderr",
    benchmark:
        "workflow/benchmarks/explode_input.txt"
    conda:
        "../envs/install.yaml"
    script:
        "../scripts/03_annotate-bases/explode.py"


# ------------- #
#     ----      #
# ------------- #


rule annotate_conservation:
    message:
        """
        Annotate control positions with conservation scores
        """
    input:
        track=rules.explode_input.output,
        linsight=LINSIGHT_SCORES,
        phastcons=PHSTCONS_SCORES,
        gerp=GERP_SCORES,
        phylop=PHYLOP_SCORES,
    output:
        Path(
            PROCESS_DIR,
            "conservation.bed",
        ),
    resources:
        mem_mb=32000,
        runtime=400,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/annotate_conservation.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/annotate_conservation.stdout",
        stderr="workflow/logs/annotate_conservation.stderr",
    script:
        "../scripts/03_annotate-bases/conservation.py"
