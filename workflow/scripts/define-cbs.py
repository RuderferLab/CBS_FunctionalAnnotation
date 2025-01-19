"""define-cbs.py"""

import pandas as pd
import pyranges as pr

# Snakemake
CCRES = snakemake.input[0]  # type:ignore
MOTIFS = snakemake.input[1]  # type:ignore
ACTIVITY = snakemake.input[2]  # type:ignore
OUTPUT = snakemake.output[0]  # type:ignore

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
        names=["Chromosome", "Start", "End", "rDHS"],
        dtype={"Chromosome": str, "Start": int, "End": int, "rDHS": str},
    )


def read_track(filepath: str) -> pd.DataFrame:
    dtype = {
        "Chromosome": str,
        "Start": int,
        "End": int,
        "score": float,
        "pval": float,
        "strand": str,
        "mid": str,
    }
    return pd.read_csv(
        filepath,
        sep="\t",
        dtype=dtype,
        engine="c",
    )


def main():
    # Read inputs
    ccres = read_ccres(CCRES)
    motifs = read_track(MOTIFS)

    # Merge activity with ccres
    activity = pd.read_csv(ACTIVITY, sep="\t", engine="c")
    ccres = ccres.merge(activity, on="rDHS", how="left")

    # Make pyRanges
    ccres = pr.PyRanges(ccres)
    motifs = pr.PyRanges(motifs)

    # Left join motifs with ccres
    ix = motifs.join(ccres).as_df()

    # Handle bookeneded overlaps by selecting the highest activity score
    ix.sort_values(["mid", "activity"], ascending=[True, False], inplace=True)
    ix.drop_duplicates("mid", keep="first", inplace=True)

    # Flag motif overlap with rDHS
    ix["is_cbs"] = 1

    # Add cbsid
    ix.reset_index(drop=True, inplace=True)
    ix["cid"] = "CBS" + (ix.index + 1).astype(str).str.zfill(6)

    # Clean up
    ix.drop(
        columns=[
            "score",
            "pval",
            "strand",
            "Start_b",
            "End_b",
        ],
        axis=1,
        inplace=True,
    )

    dtypes = {
        "Chromosome": str,
        "Start": int,
        "End": int,
        "mid": str,
        "rDHS": str,
        "activity": float,
        "quantile_rdhswide": int,
        "is_cbs": int,
        "cid": str,
    }
    ix = ix.astype(dtypes)

    # Write out
    ix.to_csv(OUTPUT, index=False, sep="\t")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
