from pathlib import Path
from snakemake.utils import min_version


# Settings
min_version("7.12.1")


# ------------- #
# Config        #
# ------------- #

# IO
INSTALL_DIR = config["IO_DIRS"]["INSTALL_DIR"]
PROCESS_DIR = config["IO_DIRS"]["PROCESS_DIR"]["STEP_5"]

# Inputs from upstream
VARIANTS = "results/04_Annotate-Variants/gnomad-snvs.processed.tsv.gz"
ATSNP = "results/04_Annotate-Variants/gnomad-snvs.atsnp.tsv.gz"
CONSERVATION = "results/03_Conservation/CBS.positions.conservation.hg38.tsv.gz"
CBSS = "results/02_Define-CBS/CBS.hg38.activity.tsv"
MOTIFS = config["IP_DATA"]["PREPROCESSED"]["CTCF_MOTIFS"]
EXPLODED_POSITIONS = "results/03_Conservation/positions.hg38.bed"

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        Path(
            PROCESS_DIR,
            "final_matrix.tsv.gz",
        ),

# ------------- #
#     ----      #
# ------------- #

rule map_cbs_positions:
    input:
        cbss=CBSS,
        exploded_posits=EXPLODED_POSITIONS,
    output:
        temp(Path(
            PROCESS_DIR,
            "cbs_positions_to_mid.tsv",
        )),
    resources:
        mem_mb=2000,
        runtime=10,
    benchmark:
        "workflow/benchmarks/map_cbs_positions.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/map_cbs_positions.stdout",
        stderr="workflow/logs/map_cbs_positions.stderr",
    shell:
        """
        vawk '{{if (NR>1) print $1, $2, $3, $4}}'  {input.cbss} |
        bedtools intersect -a {input.exploded_posits} -b stdin -loj |
        vawk '{{print $4, $8}}' > {output}
        """

# ------------- #
#     ----      #
# ------------- #


rule final_matrix:
    input:
        variants=VARIANTS,
        atsnp=ATSNP,
        conservation=CONSERVATION,
        motifs=MOTIFS,
        cbss=CBSS,
        cbs_positions_mapped=rules.map_cbs_positions.output,
    output:
        Path(
            PROCESS_DIR,
            "final_matrix.tsv.gz",
        ),
    resources:
        mem_mb=84000,
        runtime=180,
    benchmark:
        "workflow/benchmarks/final_matrix.txt"
    conda:
        "../envs/install.yaml"
    log:
        stdout="workflow/logs/final_matrix.stdout",
        stderr="workflow/logs/final_matrix.stderr",
    script:
        "../scripts/assemble-matrix.py"
