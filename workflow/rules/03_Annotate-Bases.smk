from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = Path(config["MAIN"]["IO_DIRS"]["PROCESS_DIR"]) / "03_Annotate-Bases"

# Other
ANNOTATED_MOTIFS = Path(
    "results", "02_Process-Motifs", "MA0139.2-track.masked.activity.bed"
)

# Conservation files
LINSIGHT_SCORES = config["MAIN"]["CONSERVATION"]["LINSIGHT"]
PHSTCONS_SCORES = config["MAIN"]["CONSERVATION"]["PHASTCONS"]
GERP_SCORES = config["MAIN"]["CONSERVATION"]["GERP"]
PHYLOP_SCORES = config["MAIN"]["CONSERVATION"]["PHYLOP"]

# ------------- #
# Helpers       #
# ------------- #

MOTIF_LENGTH = 15
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        Path(
            PROCESS_DIR,
            "MA0139.2-track.masked.activity.positions.bed-conservation.bed",
        ),



# ------------- #
#     ----      #
# ------------- #


rule explode_track:
    input:
        motifs=ANNOTATED_MOTIFS,
    output:
        Path(
            PROCESS_DIR,
            "MA0139.2-track.masked.activity.positions.bed",
        ),
    params:
        motif_length=MOTIF_LENGTH,
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        stdout="workflow/logs/explode_track-MA0139.2.stdout",
        stderr="workflow/logs/explode_track-MA0139.2.stderr",
    benchmark:
        "workflow/benchmarks/explode_track-MA0139.2.txt"
    conda:
        "install"
    script:
        "../scripts/03_Annotate-Bases/explode.py"


# ------------- #
#     ----      #
# ------------- #


rule annotate_conservation:
    input:
        track=rules.explode_track.output,
        linsight=LINSIGHT_SCORES,
        phastcons=PHSTCONS_SCORES,
        gerp=GERP_SCORES,
        phylop=PHYLOP_SCORES,
    output:
        Path(
            PROCESS_DIR,
            "MA0139.2-track.masked.activity.positions.bed-conservation.bed",
        ),
    resources:
        mem_mb=32000,
        runtime=400,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/annotate_conservation-MA0139.2.txt"
    conda:
        "install"
    log:
        stdout="workflow/logs/annotate_conservation-MA0139.2.stdout",
        stderr="workflow/logs/annotate_conservation-MA0139.2.stderr",
    script:
        "../scripts/03_Annotate-Bases/conservation.py"