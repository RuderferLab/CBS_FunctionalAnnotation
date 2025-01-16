from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["IO_DIRS"]["PROCESS_DIR"]["STEP_5"]

# Inputs from upstream
VARIANTS = "results/04_annotate-snvs/gnomad-snvs.processed.tsv"
ATSNP = "results/04_annotate-snvs/gnomad-snvs.atsnp.tsv"
CONSERVATION = "results/03_annotate-bases/MA0139.2-track.masked.pwm-activity.positions.conservation.bed"
POSITIONS = "results/03_annotate-bases/positions.bed"
CBSS = "results/02_process-motifs/cbs-only.noproblematic.ccre.cbsid.bed"

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        Path(
            PROCESS_DIR,
            "final_matrix.tsv.gz",
        ),


# ------------- #
#     ----      #
# ------------- #


rule final_matrix:
    input:
        variants=VARIANTS,
        atsnp=ATSNP,
        conservation=CONSERVATION,
        positions=POSITIONS,
        cbss=CBSS,
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
