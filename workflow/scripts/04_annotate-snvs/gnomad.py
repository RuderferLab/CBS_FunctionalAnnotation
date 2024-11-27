import pandas as pd
import numpy as np
from pysam import FastaFile  # type: ignore

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
SNVS = snakemake.input[1]  # type: ignore
GENOME = snakemake.input[2]  # type: ignore
MURATES = snakemake.input[3]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    """Returns variant track as df"""
    dtypes = {
            "chrm": str,
            "pos0": int,
            "pos1": int,
            "pid": str,
            "mid": str,
            "idx": int,
            "var_chrm": str,
            "var_pos0": int,
            "var_pos1": int,
            "ref": str,
            "alt": str,
            "af": str,
            "ac": str,
            "an": str,
            "vep": str,
            "caddPhred": str,
            "caddRaw": str,
            "vid": str,
    }
    return pd.read_csv(filepath, sep="\t", engine="c", header=None, names=dtypes.keys(), dtype=dtypes)

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

def mask_track(track_df: pd.DataFrame) -> pd.DataFrame:
    """Mask track"""
    return track_df.replace(".", np.nan)
    

def main():
    """d"""

    # Read track
    track = read_track(TRACK)
    
    # Mask loj
    track = mask_track(track)

    print("adding maps annotations...")

    # Setup genome
    genome_fa = genome_fasta(GENOME)

    # Add singleton flag
    track["singleton"] = (
        track["ac"].astype(float).apply(lambda x: 1 if x == 1 else 0)
    )

    # Add variant count flag if vid is not NaN
    track["isvar"] = np.where(track["vid"] == "NaN", 0, 1)

    # Annotate with tricontext
    print("Fetching trinucleotide context...")
    track["context"] = track.apply(
        lambda row: fetch_tricontext(row, genome_fa),
        axis=1,
    )

    # Explictly state NaN
    track.fillna("NaN", inplace=True)
    
    # Write output
    track.drop(columns=["var_chrm", "var_pos0", "var_pos1", "vep"], inplace=True)
    track.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
