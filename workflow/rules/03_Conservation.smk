from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["IO_DIRS"]["PROCESS_DIR"]["STEP_3"]

# CBSs defined from Step 2
CBSS = Path("results/02_Define-CBS/CBS.hg38.activity.tsv")

# HG38 to HG19 liftover chain
LIFTBACK_CHAIN = config["IP_DATA"]["UCSC"]["HG38TOHG19_CHAIN"]

# Conservation files
LINSIGHT_SCORES = config["IP_DATA"]["PREPROCESSED"]["CONSERVATION"]["LINSIGHT"]
PHSTCONS_SCORES = config["IP_DATA"]["PREPROCESSED"]["CONSERVATION"]["PHASTCONS"]
GERP_SCORES = config["IP_DATA"]["PREPROCESSED"]["CONSERVATION"]["GERP"]
PHYLOP_SCORES = config["IP_DATA"]["PREPROCESSED"]["CONSERVATION"]["PHYLOP"]


# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        Path(PROCESS_DIR, "CBS.positions.conservation.hg38.tsv"),


# ------------- #
#     ----      #
# ------------- #


rule extract_cbs_bed:
    """
    Format CBS coords into bed for downstream.
    """
    input:
        CBSS,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "cbs.bed",
            )
        ),
    resources:
        mem_mb=250,
        runtime=5,
    log:
        stdout="workflow/logs/extract_cbs_bed.stdout",
        stderr="workflow/logs/extract_cbs_bed.stderr",
    benchmark:
        "workflow/benchmarks/extract_cbs_bed.txt"
    conda:
        "../envs/install.yaml"
    shell:
        "vawk '{{if (NR>1) print $1, $2, $3}}' {input} > {output}"


# ------------- #
#     ----      #
# ------------- #


rule explode_positions:
    """
    Explodes input bed to each position in the interval.
    """
    input:
        rules.extract_cbs_bed.output,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "positions.hg38.bed",
            )
        ),
    resources:
        mem_mb=500,
        runtime=5,
    log:
        stdout="workflow/logs/explode_positions.stdout",
        stderr="workflow/logs/explode_positions.stderr",
    benchmark:
        "workflow/benchmarks/explode_positions.txt"
    conda:
        "../envs/install.yaml"
    shell:
        """
        bedops --chop 1 {input} | bedtools sort -i stdin |
        vawk '{{print $0, $1"-"$2"-"$3}}' > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule map_phastcons:
    message:
        """
        Maps exploded positions to Phastcons100 conservation scores.
        """
    input:
        track=rules.explode_positions.output,
        phastcons=PHSTCONS_SCORES,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "phastcons.hg38.bed",
            )
        ),
    resources:
        mem_mb=1000,
        runtime=15,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/annotate_conservation.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/annotate_conservation.stdout",
        stderr="workflow/logs/annotate_conservation.stderr",
    shell:
        """
        bedtools intersect -sorted -a {input.track} -b {input.phastcons} -loj |
        vawk '{{print $1, $2, $3, $4, $8}}' |
        bedtools sort -i stdin > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule map_phylop:
    message:
        """
        Maps exploded positions to PhyloP conservation scores.
        """
    input:
        track=rules.explode_positions.output,
        phylop=PHYLOP_SCORES,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "phylop.hg38.bed",
            )
        ),
    resources:
        mem_mb=1000,
        runtime=15,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/annotate_conservation.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/annotate_conservation.stdout",
        stderr="workflow/logs/annotate_conservation.stderr",
    shell:
        """
        bedtools intersect -sorted -a {input.track} -b {input.phylop} -loj |
        vawk '{{print $1, $2, $3, $4, $8}}' |
        bedtools sort -i stdin > {output}
        """

# ------------- #
#     ----      #
# ------------- #


rule liftback_positions:
    """
    Lifts back positions from hg38 to hg19. Encodes position ID to map back.
    """
    input:
        positions=rules.explode_positions.output,
        chain=LIFTBACK_CHAIN,
    output:
        mapped=temp(
            Path(
                PROCESS_DIR,
                "positions.hg19.bed",
            )
        ),
        unmapped=temp(
            Path(
                PROCESS_DIR,
                "positions.hg19.unmapped.bed",
            )
        ),
    resources:
        mem_mb=1000,
        runtime=15,
    log:
        stdout="workflow/logs/liftback_positions.stdout",
        stderr="workflow/logs/liftback_positions.stderr",
    benchmark:
        "workflow/benchmarks/liftback_positions.txt"
    conda:
        "../envs/install.yaml"
    shell:
        """
        liftOver {input.positions} {input.chain} /dev/fd/1 {output.unmapped} |
        bedtools sort -i stdin > {output.mapped}
        """


# ------------- #
#     ----      #
# ------------- #


rule map_gerp:
    message:
        """
        Maps exploded positions to GERP conservation scores.
        """
    input:
        track=rules.liftback_positions.output.mapped,
        gerp=GERP_SCORES,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "gerp.hg38.bed",
            )
        ),
    resources:
        mem_mb=1000,
        runtime=15,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/annotate_conservation.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/annotate_conservation.stdout",
        stderr="workflow/logs/annotate_conservation.stderr",
    shell:
        """
        bedtools intersect -sorted -a {input.track} -b {input.gerp} -loj |
        vawk '{{split($4, pid, "-"); print pid[1], pid[2], pid[3], $4, $8}}' |
        bedtools sort -i stdin > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule map_linsight:
    message:
        """
        Maps exploded positions to LINSIGHT conservation scores.
        """
    input:
        track=rules.liftback_positions.output.mapped,
        linsight=LINSIGHT_SCORES,
    output:
        temp(
            Path(
                PROCESS_DIR,
                "linsight.hg38.bed",
            )
        ),
    resources:
        mem_mb=1000,
        runtime=15,  # runtime in minutes
    benchmark:
        "workflow/benchmarks/annotate_conservation.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/annotate_conservation.stdout",
        stderr="workflow/logs/annotate_conservation.stderr",
    shell:
        """
        bedtools intersect -sorted -a {input.track} -b {input.linsight} -loj |
        vawk '{{split($4, pid, "-"); print pid[1], pid[2], pid[3], $4, $8}}' |
        bedtools sort -i stdin > {output}
        """


# ------------- #
#     ----      #
# ------------- #


rule annotate_conservation:
    """
    Annotates CBS positions with conservation scores.
    """
    input:
        phastcons=rules.map_phastcons.output,
        phylop=rules.map_phylop.output,
        gerp=rules.map_gerp.output,
        linsight=rules.map_linsight.output,
        positions=rules.explode_positions.output,
    output:
        Path(
            PROCESS_DIR,
            "CBS.positions.conservation.hg38.tsv",
        ),
    resources:
        mem_mb=500,
        runtime=5,
    log:
        stdout="workflow/logs/annotate_conservation.stdout",
        stderr="workflow/logs/annotate_conservation.stderr",
    benchmark:
        "workflow/benchmarks/annotate_conservation.txt"
    conda:
        "../envs/install.yaml"
    shell:
        """
        bedtools intersect -sorted -a {input.positions} -b {input.phastcons} {input.phylop} {input.gerp} {input.linsight} -loj -names phastcons phylop gerp linsight |
        vawk '{{print $4, $5, $10}}' |
        datamash -W -g1 collapse 3 |
        vawk '{{split($2, con, ","); print $1, con[1], con[2], con[3], con[4]}}' > {output}
        """
