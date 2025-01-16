"""targets.py"""

import pandas as pd
import pyranges as pr

# Snakemake
CCRES = snakemake.input[0]  # type:ignore
TRACK = snakemake.input[1]  # type:ignore
ACTIVITY = snakemake.input[2]  # type:ignore
OUTPUT_CBS = snakemake.output[0]  # type:ignore
OUTPUT_ALL = snakemake.output[1]  # type:ignore

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
    # Read CCREs and make pyranges
    ccres = read_ccres(CCRES)
    ccres = pr.PyRanges(ccres)

    # Read motif track
    track = read_track(TRACK)
    track = pr.PyRanges(track)

    # Read activity
    activity = pd.read_csv(ACTIVITY, sep="\t", engine="c")

    # Intersect track with ccres and add cluster ID to flag overlapping bookended cCREs. NOTE: left join here
    ix = track.join(ccres, how="left").cluster().as_df()

    # Merge wiht activity to help break bookended overlaps
    ix = ix.merge(activity, on="rDHS", how="left")  # NOTE: left join here

    # Handle bookeneded overlaps by selecting the highest activity score
    ix.sort_values(["Cluster", "activity"], ascending=[True, False], inplace=True)
    ix.drop_duplicates("Cluster", keep="first", inplace=True)

    # Clean up
    ix["rDHS"] = ix["rDHS"].replace("-1", ".")

    # Flag motif overlap with rDHS
    ix["CBS"] = [1 if i != "." else 0 for i in ix["rDHS"]]

    # Reduce to CBS
    cbs = ix[ix["CBS"] == 1].copy()
    cbs.reset_index(inplace=True, drop=True)

    # Add cbsid
    cbs["cid"] = "CBS" + (cbs.index + 1).astype(str).str.zfill(6)
    
    # NOTE: Adding activity quantile here
    labels = [f"{i}" for i in range(1, 101)]
    cbs["quantile_cbswide"] = pd.qcut(cbs["activity"], q=100, labels=labels)

    # Merge back with ix on index
    fields = ix.columns.to_list()
    ix = ix.merge(cbs[fields + ["cid", "quantile_cbswide"]], on=fields, how="left")

    # drop cols
    ix.drop(columns=["Start_b", "End_b", "Cluster", "CBS"], axis=1, inplace=True)

    # Clean up
    ix["cid"] = ix["cid"].fillna(".")
    dtypes = {
        "Chromosome": str,
        "Start": int,
        "End": int,
        "score": int,
        "pval": int,
        "strand": str,
        "mid": str,
        "rDHS": str,
        "activity": float,
        "quantile_rdhswide": float,
        "quantile_cbswide": float,
        "cid": str,
    }
    ix = ix.astype(dtypes)
    
    # Write out
    ix[ix['cid']!="."].to_csv(OUTPUT_CBS, index=False, sep="\t")
    ix.to_csv(OUTPUT_ALL, index=False, sep="\t")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
