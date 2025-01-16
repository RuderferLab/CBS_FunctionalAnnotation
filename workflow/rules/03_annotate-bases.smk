from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["IO_DIRS"]["PROCESS_DIR"]["STEP_3"]

# CBSs defined from Step 2
CBSS = Path("results/02_process-motifs/cbs-only.noproblematic.ccre.cbsid.bed")

# Conservation files
LINSIGHT_SCORES = config["IP_DATA"]["UCSC"]["LINSIGHT"]
PHSTCONS_SCORES = config["IP_DATA"]["UCSC"]["PHASTCONS"]
GERP_SCORES = config["IP_DATA"]["UCSC"]["GERP"]
PHYLOP_SCORES = config["IP_DATA"]["UCSC"]["PHYLOP"]


# ------------- #
# Rules         #
# ------------- #

rule all:
    input:
        Path(PROCESS_DIR, "MA0139.2-track.masked.pwm-activity.positions.conservation.bed"),


# ------------- #
#     ----      #
# ------------- #


rule explode_input:
    """
    Explodes control track into positions
    """
    input:
        CBSS,
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
            "MA0139.2-track.masked.pwm-activity.positions.conservation.bed",
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
