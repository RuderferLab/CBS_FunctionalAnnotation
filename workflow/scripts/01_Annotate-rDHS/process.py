import json
import numpy as np
import pandas as pd
from scipy import stats
from typing import Generator

# Snakemake
SIGNAL_MATRIX = snakemake.input[0]  # type: ignore
CALLED_CCRES = snakemake.input[1]  # type: ignore
SIGNAL_OUTPUT = snakemake.output[0]  # type: ignore

# Params
CHUNKSIZE = snakemake.params["chunksize"]  # type: ignore

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
        compression="gzip"
    )


def chunk_matrix(filepath: str, chunksize: int) -> Generator:
    """Returns chunked encode signal matrix"""
    for chunk in pd.read_csv(
        filepath,
        sep="\t",
        compression="gzip",
        chunksize=chunksize,
        engine="c",
        index_col="rDHS"
    ):
        yield chunk


def main():
    """Main"""

    print("Preprocessing encode signal matrix...")

    # Read CCREs
    ccres = read_ccres(CALLED_CCRES)

    # Setup storage
    dbase = []

    # Iterate over matrix in chunks
    for chunk in chunk_matrix(SIGNAL_MATRIX, CHUNKSIZE):
        
        # Tidy values
        chunk = chunk.round(4)
        
        # Calculate Stouffer's Meta Z, count is getting number of non-NaN values
        n = chunk.count(axis=1)
        meta_zscores = np.array(chunk.sum(axis=1) / np.sqrt(n))
        
        # Mask -10 values to recalc
        chunk = chunk.replace(-10, np.nan)
        
        # Calculate Stouffer's Meta Z masked
        n = chunk.count(axis=1)
        meta_zscores_masked = np.array(chunk.sum(axis=1) / np.sqrt(n))

        # Merge with CCREs into BED
        chunk = ccres.merge(chunk, how="right", on="rDHS")
        
        # Update with meta zscoress
        chunk["zscore_meta"] = meta_zscores
        chunk["zscore_meta_masked"] = meta_zscores_masked  

        # Append to dbase
        dbase.append(chunk)
       
    # Concatenate dbase
    results = pd.concat(dbase, axis=0, ignore_index=True)

    # Drop NaNs - no cooridnates
    results.dropna(subset=["chrm", "pos0", "pos1"], inplace=True)

    # Fix coord dtype
    results["pos0"] = results["pos0"].astype(int)
    results["pos1"] = results["pos1"].astype(int)
    
    # Subset
    results = results[["chrm", "pos0", "pos1", "rDHS", "zscore_meta", "zscore_meta_masked"]]

    # Save results
    results.to_csv(SIGNAL_OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
