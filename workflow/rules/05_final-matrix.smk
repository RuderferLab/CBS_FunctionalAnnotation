from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["MAIN"]["IO_DIRS"]["PROCESS_DIRS"]["STEP_5"]

# Upstream - 02_Process-Motifs
MOTIF_ACTIVITY = Path(
    "results", "02_process-motifs", "MA0139.2-track.masked.pwm-activity.bed"
)

# Upstream - 03_Annotate-Bases
TRACK_POSITIONS = Path("results", "03_annotate-bases", "positions.bed")
TRACK_CONSERVATION = Path("results", "03_annotate-bases", "conservation.bed")

# Upstream - 04_Annotate-SNVs
TRACK_VARIANTS = Path("results", "04_annotate-snvs", "gnomad-snvs.processed.tsv")
TRACK_ATSNP = Path("results", "04_annotate-snvs", "atsnp", "scored", "format", "snvs.scored.tsv")

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        Path(
            PROCESS_DIR,
            "final_matrix.tsv",
        ),


# ------------- #
#     ----      #
# ------------- #


rule final_matrix:
    input:
        track=TRACK_POSITIONS,
        activity=MOTIF_ACTIVITY,
        conservation=TRACK_CONSERVATION,
        variants=TRACK_VARIANTS,
        atsnp=TRACK_ATSNP,
    output:
        Path(
            PROCESS_DIR,
            "final_matrix.tsv",
        ),
    resources:
        mem_mb=84000,
        runtime=180,
    benchmark:
        "workflow/benchmarks/final_matrix.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/final_matrix.stdout",
        stderr="workflow/logs/final_matrix.stderr",
    script:
        "../scripts/05_final-matrix/combine.py"
