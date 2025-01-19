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
INSTALL_DIR = config["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["IO_DIRS"]["PROCESS_DIR"]["STEP_1"]

# URLs
CCRES_URL = config["URLS"]["ENCODE"]["CCRES"]
CTCF_SIGNAL_MATRIX_URL = config["URLS"]["ENCODE"]["CTCF_ZSCORE_MATRIX"]

# Test file
RDHS_SIGNAL_SUMMARY = config["IP_DATA"]["PREPROCESSED"]["RDHS_SIGNAL_SUMMARY"]

# ------------- #
#     Rules     #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "CTCF-zscore.rDHS-V3.activity.hg38.tsv"),


# ------------- #
#     ----      #
# ------------- #


rule install_ccres:
    """
    Downloads ENCODE cCRES in hg38.
    """
    output:
        Path(INSTALL_DIR, "data", "encode", "GRCh38-cCREs.bed.gz"),
    params:
        url=CCRES_URL,
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
    Downloads CTCF's ENCODE Z-score signal matrix.
    """
    output:
        Path(INSTALL_DIR, "data", "encode", "GRCh38.CTCF-zscore.rDHS-V3.txt.gz"),
    params:
        url=CTCF_SIGNAL_MATRIX_URL,
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
        ccres=rules.install_ccres.output,
        signal_summary=RDHS_SIGNAL_SUMMARY,
    output:
        Path(PROCESS_DIR, "CTCF-zscore.rDHS-V3.activity.hg38.tsv"),
    resources:
        mem_mb=2000,
        runtime=10,
    log:
        stdout=Path("workflow", "logs", "annotate_activity.stdout"),
        stderr=Path("workflow", "logs", "annotate_activity.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "annotate_activity.txt")
    conda:
        "../envs/active.yaml"
    script:
        "../scripts/rdhs-activity.py"