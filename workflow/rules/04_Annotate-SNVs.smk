from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = Path(config["MAIN"]["IO_DIRS"]["PROCESS_DIR"]) / "04_Annotate-SNVs"

# Gnomad snvs
GNOMAD_SNVS = config["MAIN"]["GNOMAD"]["SNVS"]

# JASPAR
JASPAR_PFM = Path("config", "MA0139.2-pfm.txt")

# Upstream
TRACK_POSITIONS = Path("results", "03_Annotate-Bases", "MA0139.2-track.masked.activity.positions.bed")

# Other
GENOME_DIR = Path(config["MAIN"]["GENOME"]["GENOME_DIR"])
GENOME_SIZES = GENOME_DIR / "hg38.fa.sizes"
GENOME_FASTA = GENOME_DIR / "hg38.fa"
GNOMAD_MURATES = Path("resources", "gnomad", "gnomad_v2.supplement-f10.murates-long.tsv")

# ------------ #
# Helpers       #
# ------------- #

MOTIF_LENGTH = 15
CHROMOSOMES = config["MAIN"]["TEST_PARAMS"]["CHROMOSOMES"]

       

# ------------- #
# Rules         #
# ------------- #

rule all:
    input:
        Path(
            PROCESS_DIR,
            "MA0139.2-gnomAD.v3.WGS-SNVS.bed",
        ),
        Path(
            PROCESS_DIR,
            "atSNP",
            "scored",
            "format",
            "MA0139.2_snvs.scored.tsv",
        ),


# ------------- #
#     ----      #
# ------------- #

rule annotate_snvs:
    input:
        track=TRACK_POSITIONS,
        snvs=GNOMAD_SNVS,
        genome=GENOME_FASTA,
        murates=GNOMAD_MURATES,
    output:
        Path(
            PROCESS_DIR,
            "MA0139.2-gnomAD.v3.WGS-SNVS.bed",
        ),
    log:
        stdout="workflow/logs/annotate_snvs-MA0139.2.stdout",
        stderr="workflow/logs/annotate_snvs-MA0139.2.stderr",
    resources:
        mem_mb=2000,
        runtime=60,
    benchmark:
        "workflow/benchmarks/annotate_snvs-MA0139.2.txt"
    conda:
        "install"
    script:
        "../scripts/04_Annotate-SNVs/gnomad.py"




# ------------- #
#     ----      #
# ------------- #


rule prepare_fastas:
    input:
        track=rules.annotate_snvs.output,
    output:
        ref=Path(
            PROCESS_DIR,
            "atSNP",
            "fasta",
            "ref",
            "MA0139.2-{chrom}.fasta.gz"
        ),
        alt=Path(
            PROCESS_DIR,
            "atSNP",
            "fasta",
            "alt",
            "MA0139.2-{chrom}.fasta.gz"
        ),
    params:
        slop_size=MOTIF_LENGTH - 1,
        genome_sizes=str(GENOME_SIZES),
        genome_fasta=str(GENOME_FASTA),
        chromosome=lambda wc: wc.chrom,
    log:
        stdout="workflow/logs/prepare_fastas-MA0139.2-{chrom}.stdout",
        stderr="workflow/logs/prepare_fastas-MA0139.2-{chrom}.stderr",
    resources:
        mem_mb=32000,
        runtime=60,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/prepare_fastas-MA0139.2-{chrom}.txt"
    conda:
        "install"
    script:
        "../scripts/04_Annotate-SNVs/atsnp/fasta.py"

# ------------- #
#     ----      #
# ------------- #

rule score_atsnp:
    input:
        jaspar_pfm=JASPAR_PFM,
        ref_fa=rules.prepare_fastas.output.ref,
        alt_fa=rules.prepare_fastas.output.alt,
    output:
        Path(
            PROCESS_DIR,
            "atSNP",
            "scored",
            "MA0139.2-{chrom}_snvs.scored.tsv.gz"
        ),
    log:
        stdout="workflow/logs/score_atsnp-MA0139.2-{chrom}.stdout",
        stderr="workflow/logs/score_atsnp-MA0139.2-{chrom}.stderr",
    resources:
        mem_mb=36000,
        runtime=600,
    threads: 8
    benchmark:
        "workflow/benchmarks/score_atsnp-MA0139.2-{chrom}.txt"
    conda:
        "conda_R"
    script:
        "../scripts/04_Annotate-SNVs/atsnp/score.R"

# ------------- #
#     ----      #
# ------------- #

rule format_atsnp:
    input:
        atsnp=expand(rules.score_atsnp.output, chrom=CHROMOSOMES),
    output:
        Path(
            PROCESS_DIR,
            "atSNP",
            "scored",
            "format",
            "MA0139.2_snvs.scored.tsv",
        ),
    resources:
        mem_mb=32000,
        runtime=90,
    conda:
        "install"
    log:
        stdout="workflow/logs/format_atsnp-MA0139.2.stdout",
        stderr="workflow/logs/format_atsnp-MA0139.2.stderr",
    benchmark:
        "workflow/benchmarks/format_atsnp-MA0139.2.txt"
    script:
        "../scripts/04_Annotate-SNVs/atsnp/format.py"