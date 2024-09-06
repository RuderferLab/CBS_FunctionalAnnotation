from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = Path(config["MAIN"]["IO_DIRS"]["PROCESS_DIR"]) / "05_Final-Matrix"

# Upstream - 02_Process-Motifs
MOTIF_ACTIVITY = Path(
    "results", "02_Process-Motifs", "MA0139.2-track.masked.activity.bed"
)

# Upstream - 03_Annotate-Bases
TRACK_POSITIONS = Path("results", "03_Annotate-Bases", "MA0139.2-track.masked.activity.positions.bed")
TRACK_CONSERVATION = Path(
    "results", "03_Annotate-Bases", "MA0139.2-track.masked.activity.positions.bed-conservation.bed"
)

# Upstream - 04_Annotate-SNVs
TRACK_VARIANTS = Path("results", "04_Annotate-SNVs", "MA0139.2-gnomAD.v3.WGS-SNVS.bed")
TRACK_ATSNP = Path(
    "results",
    "04_Annotate-SNVs",
    "atSNP",
    "scored",
    "format",
    "MA0139.2_snvs.scored.tsv",
)


# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        Path(
            PROCESS_DIR,
            "MA0139.2-final_matrix.tsv",
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
            "MA0139.2-final_matrix.tsv",
        ),
    resources:
        mem_mb=84000,
        runtime=180,
    benchmark:
        "workflow/benchmarks/final_matrix-MA0139.2.txt"
    conda:
        "install"
    log:
        stdout="workflow/logs/final_matrix-MA0139.2.stdout",
        stderr="workflow/logs/final_matrix-MA0139.2.stderr",
    script:
        "../scripts/05_Final_Matrix/combine.py"
