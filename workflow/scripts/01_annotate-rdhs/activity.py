import numpy as np
import pandas as pd
from typing import Generator

# Snakemake
SIGNAL_MATRIX = snakemake.input[0]  # type: ignore
CALLED_CCRES = snakemake.input[1]  # type: ignore
SIGNAL_OUTPUT = snakemake.output[0]  # type: ignore

# Params
CHUNKSIZE = 100

# ------------- #
# Functions     #
# ------------- #


def read_ccres(filepath: str) -> pd.DataFrame:
    """Returns CCRE bed file as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        engine="c",
        usecols=[0, 1, 2, 3],
        names=["chrm", "pos0", "pos1", "rDHS"],
    )


def chunk_matrix(filepath: str, chunksize: int) -> Generator:
    """Returns chunked encode signal matrix"""
    # Uncomment following if running entire dataset
    # first_chunk = True
    for chunk in pd.read_csv(
        filepath,
        sep="\t",
        index_col="rDHS",
        compression="gzip",
        chunksize=chunksize,
        engine="c",
        # skiprows=[0] if first_chunk else [],
    ):
        yield chunk
        # first_chunk = False


def calulcate_activity(rdhs_matrix: pd.DataFrame, label: str) -> pd.Series:
    """Calculate activity"""
    # Number of present (nonmissing) values
    n = rdhs_matrix.count(axis=1)
    # Stouffer's Meta Z
    return (rdhs_matrix.sum(axis=1) / np.sqrt(n)).rename(label)


def main():
    """Main"""

    print("Assigning rDHS activity...")

    # Read CCREs
    ccres = read_ccres(CALLED_CCRES)

    # Process signal matrix
    results = []
    for chunk in chunk_matrix(SIGNAL_MATRIX, CHUNKSIZE):
        # Nonmasked activity
        activity_reg = calulcate_activity(chunk, "zscore_meta")

        # Mask -10 values to recalc
        chunk = chunk.replace(-10, np.nan)

        # Masked activity
        activity_masked = calulcate_activity(chunk, "zscore_meta_masked")

        # Merge and store
        activity = pd.concat([activity_reg, activity_masked], axis=1)
        results.append(activity)

    # Concatenate dbase
    results = pd.concat(results, axis=0, ignore_index=False)

    # Merge with CCREs
    results = ccres.merge(results, how="inner", on="rDHS") # Note: inner merge for testing

    # Tidy matrix
    results = results.round(4)
    dtypes = {
        "chrm": str,
        "pos0": int,
        "pos1": int,
        "rDHS": str,
        "zscore_meta": float,
        "zscore_meta_masked": float,
    }
    results = results.astype(dtypes)

    # Save results, no header for processing later
    results.to_csv(SIGNAL_OUTPUT, sep="\t", index=False, header=None)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
