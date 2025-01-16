import pandas as pd
import numpy as np
from pysam import FastaFile  # type: ignore

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
GENOME = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    """Returns variant track as df"""
    dtypes = {
        "pid": str,
        "vid": str,
        "af": str,
        "ac": str,
        "an": str,
        "caddPhred": str,
        "caddRaw": str,
    }
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        header=None,
        names=dtypes.keys(),
        dtype=dtypes,
    )


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
    row = variant.copy()

    # Extract variant info
    if row.vid == ".":
        return "."
    else:
        info = row.vid.split("-")
        chrm = "chr" + info[0]
        pos1 = int(info[1])
        pos0 = int(pos1 - 1)

        # Set trinucleotide bounds
        l_bound = pos0 - 1
        r_bound = pos1 + 1

        # Query reference sequence for context
        return genome_sequence(chrm, l_bound, r_bound, genome)


def clean_track(track_df: pd.DataFrame) -> pd.DataFrame:
    """Mask track"""
    replace_scheme = {
        "vid": ".",
        "af": np.nan,
        "ac": np.nan,
        "an": np.nan,
        "caddPhred": np.nan,
        "caddRaw": np.nan,
    }
    dtype_scheme = {
        "pid": str,
        "vid": str,
        "af": float,
        "ac": float,
        "an": float,
        "caddPhred": float,
        "caddRaw": float,
    }
    # Replace each col
    for col, val in replace_scheme.items():
        track_df[col] = track_df[col].replace(".", val)
    return track_df.astype(dtype_scheme)


def flag_singleton(ac: float) -> float:
    if ac == 1:
        return float(1)
    elif np.isnan(ac):
        return np.nan
    else:
        return float(0)


def main():
    """d"""
    ###
    # Main
    ###

    # Read track
    track = read_track(TRACK)

    # Clean track loj
    track = clean_track(track)

    print("adding maps annotations...")

    # Setup genome
    genome_fa = genome_fasta(GENOME)

    ###
    # Add singleton flag
    ###

    track["singleton"] = track["ac"].astype(float).apply(flag_singleton)

    ###
    # Add variant count flag if vid is not "."
    ###

    track["isvar"] = np.where(track["vid"] == ".", 0, 1)

    ###
    # Annotate with tricontext
    ###

    print("Fetching trinucleotide context...")
    track["context"] = track.apply(
        lambda row: fetch_tricontext(row, genome_fa),
        axis=1,
    )

    ###
    # Write output
    ###

    track.to_csv(OUTPUT, sep="\t", index=False, compression="gzip")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
