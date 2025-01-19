"""signal-matrix.py"""

from typing import Generator
import sys
import numpy as np
import pandas as pd

# I/O
SIGNAL_MATRIX = sys.argv[1]
OUTPUT_MATRIX = "resources/files/CTCF-zscore.rDHS-V3.hg38.tsv.processed.gz"

# Parameters
CHUNKSIZE = 100000

# ------------- #
# Functions     #
# ------------- #


def chunk_matrix(filepath: str, chunksize: int) -> Generator:
    """Returns chunked encode signal matrix"""
    first_chunk = True
    for chunk in pd.read_csv(
        filepath,
        sep="\t",
        index_col="rDHS",
        compression="gzip",
        chunksize=chunksize,
        engine="c",
        skiprows=[1] if first_chunk else [],
    ):
        yield chunk


def main():
    """Main"""
    # Process signal matrix
    results = []
    for chunk in chunk_matrix(SIGNAL_MATRIX, CHUNKSIZE):

        # Mask -10 values to recalc
        chunk = chunk.replace(-10, np.nan)

        # Summarize each rDHS - sum of signal across biosamples and the number measured biosamples
        chunk["sum_signal"] = chunk.sum(axis=1)
        chunk["num_signal"] = chunk.count(axis=1)

        # Append to results
        subset_fields = ["sum_signal", "num_signal"]
        results.append(chunk[subset_fields])

    # Concatenate results
    results = pd.concat(results)

    # Extract rDHS field
    results.reset_index(inplace=True)

    # Tidy up - round to 4 decimal places and ensure dtypes
    results = results.round(4)
    results = results.astype({"rDHS": str, "sum_signal": float, "num_signal": int})

    # Write to file
    results.to_csv(OUTPUT_MATRIX, sep="\t", index=False, compression="gzip")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
