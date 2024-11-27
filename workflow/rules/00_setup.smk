from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# Signal matrix URLs
CTCF_ZSCORES_URL = config["urls"]["encode"]["signal_matrices"]["CTCF"]

# cCRE BED URL
CCRE_ALL_URL = config["urls"]["encode"]["GRCh38_cCREs"]

# Motif track and PFM URLs
JASPAR_BED_URL = config["urls"]["jaspar"]["profile_bed"]
JASPAR_PFM_URL = config["urls"]["jaspar"]["profile_pfm"]

# List of jaspar profiles to install
JASPAR_PROFILE = "MA0139.2"

# UCSC gerp++ and phastcons100 URLs
UCSC_GERP_URL = config["urls"]["ucsc"]["gerp_hg19"]
UCSC_PHYLOP100_URL = config["urls"]["ucsc"]["phylop100_hg38"]
UCSC_PHASTCONS100_URL = config["urls"]["ucsc"]["phastcons100_hg38"]

# Linsight hg19
LINSIGHT_HG19_URL = config["urls"]["others"]["linsight_hg19"]

# ------------- #
#     Rules     #
# ------------- #

rule all:
    input:
        None,




# ------------- #
#     encode    #
# ------------- #


rule encode_ccres:
    message:
        "Downloading ENCODE cCREs"
    output:
        Path("resources", "data", "encode", "GRCh38-cCREs.bed.gz"),
    params:
        url=CCRE_ALL_URL,
        download=Path("resources", "data", "encode", "GRCh38-cCREs.bed"),
    conda:
        "../envs/active.yaml"
    localrule: True
    shell:
        "curl -o {params.download} {params.url} && gzip {params.download}"


# ------------- #


rule download_encode_signal_matrix:
    message:
        "Downloading ENCODE Z-score matrix"
    output:
        Path("resources", "data", "encode", "GRCh38.CTCF-zscore.rDHS-V3.txt.gz"),
    params:
        url=lambda wc: ZSCORE_URLS[wc.matrix],
    conda:
        "../envs/active.yaml"
    localrule: True
    shell:
        "curl -o {output} {params.url}"


# ------------- #
#     jaspar    #
# ------------- #


rule jaspar_profile_pfm:
    message:
        "Downloading JASPAR profile"
    output:
        Path("resources", "data", "jaspar", "{profile}", "{profile}.jaspar"),
    params:
        profile=lambda wc: wc.profile,
        url=lambda wc: JASPAR_PFM_URL.format(profile=wc.profile),
    conda:
        "../envs/active.yaml"
    localrule: True
    # cache: "omit-software"
    shell:
        "curl -o {output} {params.url}"


# ------------- #


rule jaspar_profile_track:
    message:
        "Downloading JASPAR genome track"
    output:
        Path("resources", "data", "jaspar", "{profile}", "{profile}.tsv.gz"),
    params:
        url=lambda wc: JASPAR_PFM_URL.format(profile=wc.profile),
    conda:
        "../envs/active.yaml"
    localrule: True
    # cache: "omit-software"
    shell:
        "curl -o {output} {params.url}"


# ------------- #
#     ucsc      #
# ------------- #


rule ucsc_gerp_hg19:
    message:
        "Downloading gerp++ scores from UCSC"
    output:
        Path("resources", "data", "ucsc", "conservation", "gerp-hg19.bw"),
    params:
        url=UCSC_GERP_URL,
    conda:
        "../envs/active.yaml"
    localrule: True
    # cache: "omit-software"
    shell:
        "curl -o {output} {params.url}"


# ------------- #


rule ucsc_phastcons100_hg38:
    message:
        "Downloading phastcons100 scores from UCSC"
    output:
        Path("resources", "data", "ucsc", "conservation", "phastcons100-hg38.bw"),
    params:
        url=UCSC_PHASTCONS100_URL,
    conda:
        "../envs/active.yaml"
    localrule: True
    # cache: "omit-software"
    shell:
        "curl -o {output} {params.url}"


# ------------- #


rule ucsc_phylop100_hg38:
    message:
        "Downloading phylop100 scores from UCSC"
    output:
        Path("resources", "data", "ucsc", "conservation", "phylop100-hg38.bw"),
    params:
        url=UCSC_PHYLOP100_URL,
    conda:
        "../envs/active.yaml"
    localrule: True
    # cache: "omit-software"
    shell:
        "curl -o {output} {params.url}"


# ------------- #
#     misc      #
# ------------- #


rule linsight_hg19:
    output:
        Path("resources", "data", "other", "conservation", "linsight-hg19.bw"),
    params:
        url=LINSIGHT_HG19_URL,
    conda:
        "../envs/active.yaml"
    localrule: True
    # cache: "omit-software"
    shell:
        "curl -o {output} {params.url}"

### ** from rule 01 - should incorporate into this rule
rule install_ccres:
    """
    Downloads ENCODE ccres for human in hg38
    """
    output:
        Path(INSTALL_DIR, "data", "encode", "GRCh38-cCREs.bed"),
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
        "install"
    shell:
        "curl -o {output} {params.url}"


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
        "install"
    shell:
        "curl -o {output} {params.url}"


# ------------- #
#     ----      #
# ------------- #


rule download_profile:
    """
    Downloads JASPAR profile
    """
    output:
        Path(INSTALL_DIR, "data", "jaspar", "{profile}", "{profile}.jaspar"),
    params:
        profile=lambda wc: wc.profile,
    localrule: True
    log:
        stdout=Path("workflow", "logs", "download_profile-{profile}.stdout"),
        stderr=Path("workflow", "logs", "download_profile-{profile}.stderr"),
    conda:
        "install"
    shell:
        "curl -o {output} https://jaspar.elixir.no/api/v1/matrix/{params.profile}/?format=jaspar"


# ------------- #
#     ----      #
# ------------- #


rule download_track:
    """
    Downloads JASPAR genome track
    """
    output:
        Path(INSTALL_DIR, "data", "jaspar", "{profile}", "{profile}.tsv.gz"),
    params:
        url=JASPAR_TRACK_URL,
    localrule: True
    log:
        stdout=Path("workflow", "logs", "download_track-{profile}.stdout"),
        stderr=Path("workflow", "logs", "download_track-{profile}.stderr"),
    conda:
        "install"
    shell:
        "curl -o {output} {params.url}"