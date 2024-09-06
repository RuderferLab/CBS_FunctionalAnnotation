import pandas as pd
import numpy as np
from pysam import FastaFile  # type: ignore
import pybedtools as pbt  # type: ignore

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
SNVS = snakemake.input[1]  # type: ignore
GENOME = snakemake.input[2]  # type: ignore
MURATES = snakemake.input[3]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_murates(filepath: str) -> pd.DataFrame:
    """Returns gnomAD murate table as pandas df"""
    return pd.read_csv(filepath, sep="\t", engine="c")


def genome_fasta(genome_file: str) -> FastaFile:
    """Returns FastaFile class for input genome file"""
    return FastaFile(genome_file)


def genome_sequence(chrm: str, pos0: int, pos1: int, fasta: FastaFile) -> str:
    """Returns string of reference genome sequence for a given interval"""
    # Fetch sequence
    return fasta.fetch(chrm, pos0, pos1).upper()


def fetch_tricontext(variant: pd.Series, genome: FastaFile) -> str:
    """Returns trinucleotide context for input SNV"""
    # Work on copy and set up variant info
    info = variant.copy()
    chrm = info.chrm
    pos0 = int(info.pos0)
    pos1 = int(info.pos1)

    # Set trinucleotide bounds
    l_bound = pos0 - 1
    r_bound = pos1 + 1

    # Query reference sequence for context
    return genome_sequence(chrm, l_bound, r_bound, genome)


def mask_interval(interval):
    return ["NaN" if i == "." or i == "-1" else i for i in interval.fields]


def subset_interval(interval):
    return [interval.fields[i] for i in [0, 1, 2, 3, 4, 9, 10, 11, 12, 13, 15, 16, 17]]


def main():
    """d"""

    # Make bedtools
    track_bt = pbt.BedTool(TRACK).sort()
    snvs_bt = pbt.BedTool(SNVS)

    # Stream IX
    ix = track_bt.intersect(snvs_bt, loj=True, sorted=True, stream=True).each(
        mask_interval
    )

    print("Processing stream....")

    i = 0
    results = []
    for hit in ix:
        sub = subset_interval(hit)
        results.append(sub)

    print("Done. making df...")

    # convert to df
    results = pd.DataFrame(
        results,
        columns=[
            "chrm",
            "pos0",
            "pos1",
            "pid",
            "mid",
            "ref",
            "alt",
            "af",
            "ac",
            "an",
            "caddPhred",
            "caddRaw",
            "vid",
        ],
    )

    print("adding maps annotations...")

    # Setup genome
    genome_fa = genome_fasta(GENOME)

    # Add singleton flag
    results["singleton"] = (
        results["ac"].astype(float).apply(lambda x: 1 if x == 1 else 0)
    )

    # Add variant count flag if vid is not NaN
    results["isvar"] = np.where(results["vid"] == "NaN", 0, 1)

    # Annotate with tricontext
    print("Fetching trinucleotide context...")
    results["context"] = results.apply(
        lambda row: fetch_tricontext(row, genome_fa),
        axis=1,
    )

    # Write output
    results.drop(columns=["chrm", "pos0", "pos1"], inplace=True)
    results.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
