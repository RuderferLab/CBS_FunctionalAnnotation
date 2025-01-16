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

# Other
HG38_GAPS = config["IP_DATA"]["UCSC"]["HG38_GAPS"]
BLACKLIST = config["IP_DATA"]["ENCODE"]["BLACKLIST"]
EXONS = config["IP_DATA"]["GENCODE"]["EXONS"]
HG38_SIZES = config["IP_DATA"]["UCSC"]["HG38_SIZES"]
CCRES = config["IP_DATA"]["ENCODE"]["CCRES"]

# Activity - from section 1
ACTIVITY = Path("results/01_annotate-rdhs/GRCh38.CTCF-zscore.rDHS-V3.activity.tsv")

# ------------- #
#     Rules     #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "MA0139.2-track.masked.noproblematic.ccre.cbsid.bed"),


# ------------- #
#     ----      #
# ------------- #


rule download_profile:
    """
    Downloads JASPAR profile
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
    Downloads JASPAR genome track
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


rule build_mask:
    """
    Builds reference genome mask
    """
    input:
        gaps=HG38_GAPS,
        blacklist=BLACKLIST,
        exons=EXONS,
    output:
        Path(PROCESS_DIR, "problematic-regions.hg38.bed"),
    resources:
        mem_mb=2000,
        runtime=30,
    log:
        stdout=Path("workflow", "logs", "build_mask"),
        stderr=Path("workflow", "logs", "build_mask"),
    benchmark:
        Path("workflow", "benchmarks", "build_mask")
    conda:
        "../envs/install.yaml"
    shell:
        """
        zcat {input.gaps} {input.blacklist} {input.exons} | 
        vawk '{{print $1, $2, $3}}' | 
        sort -k1,1 -k2,2n |
        bedtools merge -i stdin > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule mask_track:
    """
    Removes excluded regions from motif track
    """
    input:
        track=rules.download_track.output,
        mask=rules.build_mask.output,
    output:
        temp(Path(PROCESS_DIR, "MA0139.2-track.masked.tsv")),
    resources:
        mem_mb=2000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "mask_track-MA0139.2.stdout"),
        stderr=Path("workflow", "logs", "mask_track-MA0139.2.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "jaspar", "mask_track-MA0139.2.txt")
    conda:
        "../envs/install.yaml"
    shell:
        "bedtools intersect -a {input.track} -b {input.mask} -v > {output}"


# ------------- #
#     ----      #
# ------------- #


rule remove_track_problematic:
    """
    Removes problematic motifs from track
    """
    input:
        track=rules.mask_track.output,
    output:
        Path(PROCESS_DIR, "MA0139.2-track.masked.noproblematic.bed"),
    resources:
        mem_mb=2000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "mask_track-MA0139.2.stdout"),
        stderr=Path("workflow", "logs", "mask_track-MA0139.2.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "jaspar", "mask_track-MA0139.2.txt")
    conda:
        "../envs/install.yaml"
    script:
        "../scripts/02_process-motifs/problematic.py"


# ------------- #
#     ----      #
# ------------- #

rule annotate_cbs:
    """
    Annotate motifs with CBS ID
    """
    input:
        ccres=CCRES,
        track=rules.remove_track_problematic.output,
        activity=ACTIVITY,
    output:
        cbs=Path(PROCESS_DIR, "cbs-only.noproblematic.ccre.cbsid.bed"),
        compete=Path(PROCESS_DIR, "MA0139.2-track.masked.noproblematic.ccre.cbsid.bed"),
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
       "../scripts/02_process-motifs/targets.py"