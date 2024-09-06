from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = Path(config["MAIN"]["IO_DIRS"]["PROCESS_DIR"]) / "02_Process-Motifs"

# MOTIF TRACK URL
JASPAR_TRACK = config["MAIN"]["TEST_FILES"]["JASPAR_TRACK"]

# Genome
GENOME_DIR = Path(config["MAIN"]["GENOME"]["GENOME_DIR"])
GENOME_GAPS = GENOME_DIR / "hg38.gaps.bed"
BLACKLIST = GENOME_DIR / "hg38.blacklist.bed"
EXONS = Path(config["MAIN"]["GENOME"]["EXONS"])

# Activity - from section 1
ACTIVITY = Path("results/01_Annotate-rDHS/GRCh38-cCRES.CTCF_signal.process.tsv")

# ------------- #
#     Rules     #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "MA0139.2-track.masked.activity.bed"),


# ------------- #
#     ----      #
# ------------- #


rule mask_track:
    input:
        track=JASPAR_TRACK,
        gaps=GENOME_GAPS,
        blacklist=BLACKLIST,
        exons=EXONS,
    output:
        temp(Path(PROCESS_DIR, "MA0139.2-track.masked.bed.gz")),
    resources:
        mem_mb=2000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "mask_track-MA0139.2.stdout"),
        stderr=Path("workflow", "logs", "mask_track-MA0139.2.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "jaspar", "mask_track-MA0139.2.txt")
    conda:
        "install"
    script:
        "../scripts/02_Process-Motifs/problematic.py"


# ------------- #
#     ----      #
# ------------- #


rule annotate_activity:
    """
    Annotated motifs with activity scores from section 1
    """
    input:
        motifs=rules.mask_track.output,
        activity=ACTIVITY,
    output:
        matrix=Path(PROCESS_DIR, "MA0139.2-track.masked.activity.bed"),
    resources:
        mem_mb=2000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "annotate_activity-MA0139.2.stdout"),
        stderr=Path("workflow", "logs", "annotate_activity-MA0139.2.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "annotate_activity-MA0139.2.txt")
    conda:
        "install"
    script:
        "../scripts/02_Process-Motifs/annotate.py"

