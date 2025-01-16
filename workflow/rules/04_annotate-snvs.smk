from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["IO_DIRS"]["PROCESS_DIR"]["STEP_4"]

# Gnomad snvs
GNOMAD_SNVS = config["IP_DATA"]["GNOMAD"]["SNVS"]

# JASPAR
JASPAR_PFM = Path("resources", "data", "jaspar", "MA0139.2-pfm.txt")

# Upstream
TRACK_POSITIONS = Path("results", "03_annotate-bases", "positions.bed")

# Other
HG38_SIZES = config["IP_DATA"]["UCSC"]["HG38_SIZES"]
HG38_FASTA = config["IP_DATA"]["UCSC"]["HG38_FASTA"]
GNOMAD_MURATES = Path(
    "resources", "data", "gnomad", "gnomad_v2.supplement-f10.murates-long.tsv"
)


# ------------- #
# Rules         #
# ------------- #

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

rule all:
    input:
        Path(PROCESS_DIR, "gnomad-snvs.atsnp.tsv"),
        Path(PROCESS_DIR, "gnomad-snvs.processed.tsv")


# ------------- #
#     ----      #
# ------------- #


rule format_positions:
    """
    Formats positinos for intesect
    """
    input:
        TRACK_POSITIONS,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "positions.for_ix.bed",
            )
        ),
    log:
        stdout="workflow/logs/format_positions.stdout",
        stderr="workflow/logs/format_positions.stderr",
    resources:
        mem_mb=16000,
        runtime=120,
    benchmark:
        "workflow/benchmarks/format_positions.txt"
    conda:
        "../envs/install.yaml"
    shell:
        """
        vawk '{{split($1,pid,"-"); print pid[1], pid[2]-1, pid[2], $1}}' {input} | bedtools sort -i stdin > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule extract_snvs:
    """
    Extract SNVs from gnomAD
    """
    input:
        track=rules.format_positions.output,
        snvs=GNOMAD_SNVS,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "gnomad-snvs.tsv",
            )
        ),
    log:
        stdout="workflow/logs/extract_snvs.stdout",
        stderr="workflow/logs/extract_snvs.stderr",
    resources:
        mem_mb=16000,
        runtime=120,
    benchmark:
        "workflow/benchmarks/extract_snvs.txt"
    conda:
        "../envs/install.yaml"
    shell:
        """
        bedtools intersect -a {input.track} -b {input.snvs} -loj -sorted |
        vawk '{{print $4, $16, $10, $11, $12, $14, $15}}' > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule process_snvs:
    """
    Annotates SNVs with maps annotations
    """
    input:
        track=rules.extract_snvs.output,
        genome=HG38_FASTA,
    output:
        Path(
            PROCESS_DIR,
            "gnomad-snvs.processed.tsv.gz",
        ),
    log:
        stdout="workflow/logs/process_snvs.stdout",
        stderr="workflow/logs/process_snvs.stderr",
    resources:
        mem_mb=16000,
        runtime=60,
    benchmark:
        "workflow/benchmarks/process_snvs.txt"
    conda:
        "../envs/install.yaml"
    script:
        "../scripts/04_annotate-snvs/gnomad.py"


# ------------- #
#     ----      #
# ------------- #

rule prepare_fastas:
    """
    Prepare reference and alternative fasta files for atSNP
    """
    input:
        track=rules.process_snvs.output,
    output:
        ref=expand(Path(PROCESS_DIR, "atsnp", "fasta", "ref", "{chrom}.fasta"), chrom=CHROMOSOMES),
        alt=expand(Path(PROCESS_DIR, "atsnp", "fasta", "alt", "{chrom}.fasta"), chrom=CHROMOSOMES),
    params:
        slop_size=14,
        genome_sizes=str(HG38_SIZES),
        genome_fasta=str(HG38_FASTA),
        outdir=f"{PROCESS_DIR}/atsnp/fasta",
    log:
        stdout="workflow/logs/prepare_fastas.stdout",
        stderr="workflow/logs/prepare_fastas.stderr",
    resources:
        mem_mb=32000,
        runtime=90,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/prepare_fastas.txt"
    conda:
        "../envs/install.yaml"
    script:
        "../scripts/04_annotate-snvs/atsnp/fasta.py"


# ------------- #
#     ----      #
# ------------- #


rule score_atsnp:
    """
    Implement atSNP scoring
    """
    input:
        jaspar_pfm=JASPAR_PFM,
        ref_fa=Path(PROCESS_DIR, "atsnp", "fasta", "ref", "{chrom}.fasta.gz"),
        alt_fa=Path(PROCESS_DIR, "atsnp", "fasta", "alt", "{chrom}.fasta.gz"),
    output:
        Path(
            PROCESS_DIR,
            "atsnp",
            "scored",
            "chroms",
            "{chrom}-snvs.scored.tsv.gz",
        ),
    log:
        stdout="workflow/logs/score_atsnp-{chrom}.stdout",
        stderr="workflow/logs/score_atsnp-{chrom}.stderr",
    resources:
        mem_mb=96000,
        runtime=600,  # runtime in minutes
    threads: 8
    benchmark:
        "workflow/benchmarks/score_atsnp-{chrom}.txt"
    conda:
        "conda_R" # TODO: add/troubleshoot conda env
    script:
        "../scripts/04_annotate-snvs/atsnp/score.R"


# ------------- #
#     ----      #
# ------------- #

rule format_atsnp:
    message:
        """
        Combine and format atSNP scores
        """
    input:
        track=rules.process_snvs.output,
        atsnp=expand(rules.score_atsnp.output, chrom=CHROMOSOMES),
    output:
        Path(
            PROCESS_DIR,
            "gnomad-snvs.atsnp.tsv.gz",
        ),
    resources:
        mem_mb=32000,
        runtime=90,  # runtime in minutes
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/format_atsnp.stdout",
        stderr="workflow/logs/format_atsnp.stderr",
    benchmark:
        "workflow/benchmarks/format_atsnp-.txt"
    script:
        "../scripts/04_annotate-snvs/atsnp/format.py"
