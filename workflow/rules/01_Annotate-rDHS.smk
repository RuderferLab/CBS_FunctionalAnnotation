from pathlib import Path
from snakemake.utils import min_version

# ------------- #
# Settings
# ------------- #

min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = Path(config["MAIN"]["IO_DIRS"]["PROCESS_DIR"]) / "01_Annotate-rDHS"

# URLs
ENCODE_CCRES = config["MAIN"]["TEST_FILES"]["ENCODE_CCRES"]
CTCF_SIGNAL_MATRIX_SAMPLE = config["MAIN"]["TEST_FILES"]["CTCF_ZSCORE_MATRIX_SAMPLE"]

# ------------- #
#     Rules     #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "GRCh38-cCRES.CTCF_signal.process.tsv"),
        #Path(PROCESS_DIR, "GRCh38-cCRES.CTCF_signal.activity.tsv"),

# ------------- #
#     ----      #
# ------------- #


rule process_signal_matrx:
    input:
        signal_matrix=CTCF_SIGNAL_MATRIX_SAMPLE,
        called_ccres=ENCODE_CCRES,
    output:
        Path(PROCESS_DIR, "GRCh38-cCRES.CTCF_signal.process.tsv"),
    params:
        chunksize=1000,
    resources:
        mem_mb=12000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "process_signal_matrx.stdout"),
        stderr=Path("workflow", "logs", "process_signal_matrx.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "process_signal_matrx.txt")
    conda:
        "install"
    script:
        "../scripts/01_Annotate-rDHS/process.py"