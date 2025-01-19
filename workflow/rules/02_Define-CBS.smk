import pandas as pd
from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["IO_DIRS"]["PROCESS_DIR"]["STEP_2"]

# MOTIF TRACK URL
JASPAR_TRACK_URL = config["URLS"]["JASPAR"]["TRACK"]

# ENCODE cCREs
CCRES = config["IP_DATA"]["ENCODE"]["CCRES"]

# Preprocessed motifs
MOTIFS = config["IP_DATA"]["PREPROCESSED"]["CTCF_MOTIFS"]

# Activity - from section 1
ACTIVITY = Path("results/01_rDHS-Activity/CTCF-zscore.rDHS-V3.activity.hg38.tsv")

# ------------- #
#     Rules     #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "CBS.hg38.activity.tsv"),


# ------------- #
#     ----      #
# ------------- #


rule download_profile:
    """
    Downloads CTCF's JASPAR profile MA0139.2.
    """
    output:
        Path(INSTALL_DIR, "data", "jaspar", "MA0139.2.jaspar"),
    localrule: True
    log:
        stdout=Path("workflow", "logs", "download_profile.stdout"),
        stderr=Path("workflow", "logs", "download_profile.stderr"),
    conda:
        "../envs/install.yaml"
    shell:
        "curl -o {output} https://jaspar.elixir.no/api/v1/matrix/MA0139.2/?format=jaspar"


# ------------- #
#     ----      #
# ------------- #


rule download_track:
    """
    Downloads CTCF's PWM genome track from JASPAR.
    """
    output:
        Path(INSTALL_DIR, "data", "jaspar", "MA0139.2.tsv.gz"),
    params:
        url=JASPAR_TRACK_URL,
    localrule: True
    log:
        stdout=Path("workflow", "logs", "download_track-MA0139.2.stdout"),
        stderr=Path("workflow", "logs", "download_track-MA0139.2.stderr"),
    conda:
        "../envs/install.yaml"
    shell:
        "curl -o {output} {params.url}"


# ------------- #
#     ----      #
# ------------- #

rule define_cbs:
    """
    Defines CBSs using ENCODE and JASPAR datasets; annotates activity.
    """
    input:
        ccres=CCRES,
        track=MOTIFS,
        activity=ACTIVITY,
    output:
        Path(PROCESS_DIR, "CBS.hg38.activity.tsv"),
    resources:
        mem_mb=2000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "annotate_cbs.stdout"),
        stderr=Path("workflow", "logs", "annotate_cbs.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "annotate_cbs.txt")
    conda:
        # "../envs/install.yaml"
        "install"
    script:
       "../scripts/define-cbs.py"