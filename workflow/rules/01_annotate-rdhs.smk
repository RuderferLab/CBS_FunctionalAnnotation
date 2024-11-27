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
PROCESS_DIR = config["MAIN"]["IO_DIRS"]["PROCESS_DIRS"]["STEP_1"]

# URLs
CCRE_CALLS = config["MAIN"]["URLS"]["ENCODE"]["CCRES"]
CTCF_SIGNAL_MATRIX = config["MAIN"]["URLS"]["ENCODE"]["ZSCORES"]

# Test file
CTCF_SIGNAL_MATRIX_TEST = config["MAIN"]["TEST_FILES"]["CTCF_ZSCORE_MATRIX_SAMPLE"]

# ------------- #
#     Rules     #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "rdhs-activity.hg38.tsv"),


# ------------- #
#     ----      #
# ------------- #


rule install_ccres:
    """
    Downloads ENCODE ccres for human in hg38
    """
    output:
        Path(INSTALL_DIR, "data", "encode", "GRCh38-cCREs.bed.gz"),
    params:
        url=CCRE_CALLS,
    resources:
        mem_mb=150,
        runtime=30,
    log:
        stdout=Path("workflow", "logs", "encode", "install_ccres.stdout"),
        stderr=Path("workflow", "logs", "encode", "install_ccres.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "encode", "install_ccres.txt")
    conda:
        "../envs/install.yaml"
    shell:
        "curl {params.url} | gzip > {output}"


# ------------- #
#     ----      #
# ------------- #


rule install_signal_matrix:
    """
    Downloads ENCODE Z-score matrices
    """
    output:
        Path(INSTALL_DIR, "data", "encode", "GRCh38.CTCF-zscore.rDHS-V3.txt.gz"),
    params:
        url=CTCF_SIGNAL_MATRIX,
    resources:
        mem_mb=250,
        runtime=75,
    log:
        stdout=Path("workflow", "logs", "encode", "install_signal_matrix.stdout"),
        stderr=Path("workflow", "logs", "encode", "install_signal_matrix.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "encode", "install_signal_matrix.txt")
    conda:
        "../envs/install.yaml"
    shell:
        "curl -o {output} {params.url}"


# ------------- #
#     ----      #
# ------------- #


rule annotate_activity:
    """
    Annotates rDHS with activity scores.
    """
    input:
        signal_matrix=CTCF_SIGNAL_MATRIX_TEST,
        called_ccres=rules.install_ccres.output,
    output:
        Path(PROCESS_DIR, "rdhs-activity.hg38.tsv"),
    params:
        chunksize=1000,
    resources:
        mem_mb=12000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "annotate_activity.stdout"),
        stderr=Path("workflow", "logs", "annotate_activity.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "annotate_activity.txt")
    conda:
        "../envs/active.yaml"
    script:
        "../scripts/01_annotate-rdhs/activity.py"