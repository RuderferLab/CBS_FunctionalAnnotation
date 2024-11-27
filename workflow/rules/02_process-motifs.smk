from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["MAIN"]["IO_DIRS"]["PROCESS_DIRS"]["STEP_2"]

# MOTIF TRACK URL
JASPAR_TRACK_URL = config["MAIN"]["TEST_FILES"]["JASPAR_TRACK"]

# Other
GENOME_GAPS = config["MAIN"]["GENOME"]["GENOME_GAPS"]
BLACKLIST = config["MAIN"]["GENOME"]["BLACKLIST"]
EXONS = config["MAIN"]["GENOME"]["EXONS"]
CHROM_SIZES = config["MAIN"]["GENOME"]["CHROM_SIZES"]
CCRES = config["MAIN"]["TEST_FILES"]["ENCODE_CCRES"]

# Activity - from section 1
ACTIVITY = Path("results/01_annotate-rdhs/rdhs-activity.hg38.tsv")

# ------------- #
#     Rules     #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "hg38.masked_regions.bed"),
        Path(PROCESS_DIR, "MA0139.2-track.masked.bed"),
        Path(PROCESS_DIR, "MA0139.2-track.masked.pwm-activity.bed"),


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
        gaps=GENOME_GAPS,
        blacklist=BLACKLIST,
        exons=EXONS,
    output:
        Path(PROCESS_DIR, "hg38.masked_regions.bed"),
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
        "zcat {input.gaps} {input.blacklist} {input.exons} | vawk '{{print $1, $2, $3}}' | sort -k1,1 -k2,2n > {output}"


# ------------- #
#     ----      #
# ------------- #


rule mask_track:
    """
    Masks JASPAR track based on reference genome mask
    """
    input:
        track=rules.download_track.output,
        mask=rules.build_mask.output,
    output:
        Path(PROCESS_DIR, "MA0139.2-track.masked.bed"),
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


rule make_targets:
    """
    Annotated motifs with activity scores from section 1
    """
    input:
        motifs=rules.mask_track.output,
        activity=ACTIVITY,
    output:
        Path(PROCESS_DIR, "MA0139.2-track.masked.pwm-activity.bed"),
    resources:
        mem_mb=2000,
        runtime=45,
    log:
        stdout=Path("workflow", "logs", "make_targets.stdout"),
        stderr=Path("workflow", "logs", "make_targets.stderr"),
    benchmark:
        Path("workflow", "benchmarks", "make_targets.txt")
    conda:
        "../envs/install.yaml"
    shell:
        """
        bedtools intersect -a {input.activity} -b {input.motifs} | 
        bedtools intersect -a stdin -b {input.motifs} -loj |
        vawk '{{print $7, $8, $9, $10, $11, $12, $13, $14, $4, $5, $6}}' > {output}
        """
