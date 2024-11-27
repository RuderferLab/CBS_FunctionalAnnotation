from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["MAIN"]["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["MAIN"]["IO_DIRS"]["PROCESS_DIRS"]["STEP_4"]

# Gnomad snvs
GNOMAD_SNVS = config["MAIN"]["GNOMAD"]["SNVS"]

# JASPAR
JASPAR_PFM = Path("resources", "data", "jaspar", "MA0139.2-pfm.txt")

# Upstream
TRACK_POSITIONS = Path("results", "03_annotate-bases", "positions.bed")

# Other
GENOME_SIZES = config["MAIN"]["GENOME"]["CHROM_SIZES"]
GENOME_FASTA =  config["MAIN"]["GENOME"]["HG38_FASTA"]
GNOMAD_MURATES = Path("resources", "data", "gnomad", "gnomad_v2.supplement-f10.murates-long.tsv")


# ------------- #
# Rules         #
# ------------- #

# CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"] NOTE: for full dataset
CHROMOSOMES = ["chr1", "chr10", "chr11", "chr12"]

rule all:
    input:
        Path(
            PROCESS_DIR,
            "atsnp",
            "scored",
            "format",
            "snvs.scored.tsv",
        ),



# ------------- #
#     ----      #
# ------------- #


rule extract_snvs:
    """
    Extract SNVs from gnomAD
    """
    input:
        track=TRACK_POSITIONS,
        snvs=GNOMAD_SNVS,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "gnomad-snvs.tsv",
            )
        ),
    log:
        stdout="workflow/logs/annotate_snvs.stdout",
        stderr="workflow/logs/annotate_snvs.stderr",
    resources:
        mem_mb=16000,
        runtime=120,
    benchmark:
        "workflow/benchmarks/annotate_snvs.txt"
    conda:
        "../envs/install.yaml"
    shell:
        """
        bedtools sort -i {input.track} |
        bedtools intersect -a stdin -b {input.snvs} -loj -sorted > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule process_snvs:
    """
    Annotates SNVs with murate, tricontext
    """
    input:
        track=rules.extract_snvs.output,
        snvs=GNOMAD_SNVS,
        genome=GENOME_FASTA,
        murates=GNOMAD_MURATES,
    output:
        Path(
            PROCESS_DIR,
            "gnomad-snvs.processed.tsv",
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
        ref=Path(PROCESS_DIR, "atsnp", "fasta", "ref", "{chrom}.fasta.gz"),
        alt=Path(PROCESS_DIR, "atsnp", "fasta", "alt", "{chrom}.fasta.gz"),
    params:
        slop_size=14,
        genome_sizes=str(GENOME_SIZES),
        genome_fasta=str(GENOME_FASTA),
        chromosome=lambda wc: wc.chrom,
    log:
        stdout="workflow/logs/prepare_fastas-{chrom}.stdout",
        stderr="workflow/logs/prepare_fastas-{chrom}.stderr",
    resources:
        mem_mb=32000,
        runtime=90,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/prepare_fastas-{chrom}.txt"
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
        ref_fa=rules.prepare_fastas.output.ref,
        alt_fa=rules.prepare_fastas.output.alt,
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
        atsnp=expand(rules.score_atsnp.output, chrom=CHROMOSOMES),
    output:
        Path(
            PROCESS_DIR,
            "atsnp",
            "scored",
            "format",
            "snvs.scored.tsv",
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
